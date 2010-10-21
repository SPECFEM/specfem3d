
  program read_gocad_block_extract

! when compiling this code on a PC, use "pgf90 -byteswapio " to byte-swap
! the original binary file, which was created on an SGI

  implicit none

  include "../../../constants.h"
  include "../../../constants_gocad.h"

!
! new hi-res Voxet Peter July 29, 2002
!


! be careful: fictitious velocity is different in this high-res block
!  (-99999. km/s  instead of +20. km/s)


! AXIS_O 371052.25 3774000 400
! AXIS_U 46000 0 0
! AXIS_V 0 -48750 0
! AXIS_W 0 0 -9900
! AXIS_MIN 0 0 0
! AXIS_MAX 1 1 1
! AXIS_N 185 196 100
!

  real vp_block_gocad_sngl(0:NX_GOCAD_HR-1,0:NY_GOCAD_HR-1,0:NZ_GOCAD_HR-1)
  double precision vp_block_gocad(0:NX_GOCAD_HR-1,0:NY_GOCAD_HR-1,0:NZ_GOCAD_HR-1)
  logical iflag_point(0:NX_GOCAD_HR-1,0:NY_GOCAD_HR-1,0:NZ_GOCAD_HR-1)
  integer ipoin_store(0:NX_GOCAD_HR-1,0:NY_GOCAD_HR-1,0:NZ_GOCAD_HR-1)

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
  open(unit=27,file='LA_HR_VINT_clean@@',status='old',access='direct',recl=4)
  irecord = 0
  do iz = 0,NZ_GOCAD_HR-1
  do iy = 0,NY_GOCAD_HR-1
  do ix = 0,NX_GOCAD_HR-1
    irecord = irecord + 1
    read(27,rec=irecord) vp_block_gocad_sngl(ix,iy,iz)
  enddo
  enddo
  enddo
  close(27)

! convert block read to double precision
  vp_block_gocad(:,:,:) = dble(vp_block_gocad_sngl(:,:,:))

! extend basin model below threshold to bottom of the grid to make sure
! there is no small gap between interpolated basement map and sediments
  if(EXTEND_VOXET_BELOW_BASEMENT) then
  do ix = 0,NX_GOCAD_HR-1
    do iy = 0,NY_GOCAD_HR-1

! determine at what depth the sediments, if any, end
! a fictitious P velocity of -9999 km/s has been used to flag fictitious points
      iz_found = -1
      do iz = NZ_GOCAD_HR-1,0,-1
        if(vp_block_gocad(ix,iy,iz) < 6499.) iz_found = iz
      enddo

! if some sediments are detected on this vertical line in Voxet
      if(iz_found > -1) then

! define Gocad grid, shift of Voxet is taken into account
        zsedim_found = ORIG_Z_GOCAD_HR + iz_found*SPACING_Z_GOCAD_HR

! if point is below threshold, we know we honor it with the mesh,
! therefore we can safely extend below to make sure we leave no small gap
! between our mesh and the Gocad voxet (because we interpolate the basement
! slightly differently from what has been done in Gocad at Harvard)
        if(zsedim_found <= Z_THRESHOLD_HONOR_BASEMENT) then
          do iz = max(1,iz_found-NCELLS_EXTEND),iz_found-1
            vp_block_gocad(ix,iy,iz) = vp_block_gocad(ix,iy,iz_found)
          enddo
        endif

      endif

    enddo
  enddo
  endif

!---

! also make sure there are no gaps between topography and sediments
! because we also define topography slightly differently from Gocad

  if(EXTEND_VOXET_ABOVE_TOPO) then

  print *,'reading topography from file to fill small gaps'
  call read_basin_topo_bathy_file(itopo_bathy_basin)
  print *,'regional topography file read ranges in m from ', &
             minval(itopo_bathy_basin),' to ',maxval(itopo_bathy_basin)
  print *

  max_error = -1000000.

  do ix = 0,NX_GOCAD_HR-1
    do iy = 0,NY_GOCAD_HR-1

! determine at what height the sediments, if any, end
! a fictitious P velocity of -9999 km/s has been used to flag fictitious points
      iz_found = -1
      do iz = 0,NZ_GOCAD_HR-1
        if(vp_block_gocad(ix,iy,iz) < 6499.) iz_found = iz
      enddo

! if some sediments are detected on this vertical line in Voxet
      if(iz_found > -1) then

! define Gocad grid, shift of Voxet is taken into account
        xcoord = ORIG_X_GOCAD_HR + ix*SPACING_X_GOCAD_HR
        ycoord = ORIG_Y_GOCAD_HR + iy*SPACING_Y_GOCAD_HR
        zsedim_found = ORIG_Z_GOCAD_HR + iz_found*SPACING_Z_GOCAD_HR

! determine height of topography at this point
  call utmgeo(long,lat,xcoord,ycoord,IZONE_UTM_LA,IUTM2LONGLAT)

! get closest point in bathy/topo model
  iclosestlong = nint((long - ORIG_LONG) / DEGREES_PER_CELL) + 1
  iclosestlat = nint((lat - ORIG_LAT) / DEGREES_PER_CELL) + 1

! avoid edge effects and extend with identical topo if point outside model
  if(iclosestlong < 1) iclosestlong = 1
  if(iclosestlong > NX_TOPO) iclosestlong = NX_TOPO

  if(iclosestlat < 1) iclosestlat = 1
  if(iclosestlat > NY_TOPO) iclosestlat = NY_TOPO

! compute elevation at current point
    elevation = dble(itopo_bathy_basin(iclosestlong,iclosestlat))

! if distance is negative, it means our topo is below Gocad topo
! compute maximum to estimate maximum error between the two surfaces
    if(elevation - zsedim_found < 0.d0) max_error = dmax1(max_error,dabs(elevation - zsedim_found))

! if point is not too far from topo, assume sediments should reach the surface,
! and fill the gap and extend above topo to be safe
    if(elevation - zsedim_found < DISTMAX_ASSUME_SEDIMENTS) then
      do iz = iz_found+1,min(iz_found+NCELLS_EXTEND,NZ_GOCAD_HR-1)
        vp_block_gocad(ix,iy,iz) = vp_block_gocad(ix,iy,iz_found)
      enddo
    endif

      endif

    enddo
  enddo

  print *
  print *,'maximum error detected between our topo and Gocad = ',max_error
  print *

  endif


!---

  icount_undefined = 0
  vpmin = + 100000000.
  vpmax = - 100000000.
  ipoin = 0

! count total number of points kept
  do ix = 0,NX_GOCAD_HR-1
    do iy = 0,NY_GOCAD_HR-1
      do iz = 0,NZ_GOCAD_HR-1

! exclude points that are undefined
! a fictitious P velocity of 6501 km/s has been used to flag these points
!!!! DK DK UGLY CRADE   ugly to extract only one layer for AVS
!!!! DK DK UGLY CRADE    if(vp_block_gocad(ix,iy,iz) > 6499. .or. (iz /= 90 .and. iz /= 91)) then

        if(vp_block_gocad(ix,iy,iz) > 6499.) then
          icount_undefined = icount_undefined + 1
          iflag_point(ix,iy,iz) = .false.

        else
          ipoin = ipoin + 1
          ipoin_store(ix,iy,iz) = ipoin
          iflag_point(ix,iy,iz) = .true.

          vpmin = dmin1(vpmin,vp_block_gocad(ix,iy,iz))
          vpmax = dmax1(vpmax,vp_block_gocad(ix,iy,iz))

        endif

      enddo
    enddo
  enddo

  npoin = ipoin

! export model to text file for portability
  open(unit=27,file='LA_HR_voxet_extracted.txt',status='unknown')

! write total number of points
  write(27,*) npoin

! write point found to text file, ignore decimals for Vp
! exclude points that are undefined
  do ix = 0,NX_GOCAD_HR-1
    do iy = 0,NY_GOCAD_HR-1
      do iz = 0,NZ_GOCAD_HR-1
        if(iflag_point(ix,iy,iz)) write(27,*) ix,' ',iy,' ',iz,' ',nint(vp_block_gocad(ix,iy,iz))
      enddo
    enddo
  enddo

  close(27)

! count total number of elements
  ispec = 0
  do ix = 0,NX_GOCAD_HR-2
    do iy = 0,NY_GOCAD_HR-2
      do iz = 0,NZ_GOCAD_HR-2

! suppress elements that are undefined
   if(iflag_point(ix,iy,iz) .and. &
      iflag_point(ix+1,iy,iz) .and. &
      iflag_point(ix+1,iy+1,iz) .and. &
      iflag_point(ix,iy+1,iz) .and. &
      iflag_point(ix,iy,iz+1) .and. &
      iflag_point(ix+1,iy,iz+1) .and. &
      iflag_point(ix+1,iy+1,iz+1) .and. &
      iflag_point(ix,iy+1,iz+1)) ispec = ispec + 1

      enddo
    enddo
  enddo

  nspec = ispec

  print *
  print *,'found ',icount_undefined,' undefined points'
  print *,'which is ',100.*icount_undefined/dble(NX_GOCAD_HR*NY_GOCAD_HR*NZ_GOCAD_HR),' %'
  print *

  print *,'minval maxval Vp in region kept = ',vpmin,vpmax

! create AVS file with velocity block
! only save elements that are not fictitious

  print *
  print *,'creating AVS file'
  print *,'points = ',npoin
  print *,'elements = ',nspec
  print *

  open(unit=11,file='AVS_Pvelocity_block_HR.inp',status='unknown')

! write AVS header with element data
  write(11,*) npoin,' ',nspec,' 1 0 0'

! output global AVS points
! loop on all the points
  ipoin = 0
  do ix = 0,NX_GOCAD_HR-1
    do iy = 0,NY_GOCAD_HR-1
      do iz = 0,NZ_GOCAD_HR-1

    if(iflag_point(ix,iy,iz)) then
        ipoin = ipoin + 1

! define Gocad grid, shift of Voxet is taken into account
        xcoord = ORIG_X_GOCAD_HR + ix*SPACING_X_GOCAD_HR
        ycoord = ORIG_Y_GOCAD_HR + iy*SPACING_Y_GOCAD_HR
        zcoord = ORIG_Z_GOCAD_HR + iz*SPACING_Z_GOCAD_HR
        write(11,*) ipoin,' ',sngl(xcoord),' ',sngl(ycoord),' ',sngl(zcoord)
     endif

      enddo
    enddo
  enddo

! output global AVS elements
! loop on all the elements
  ispec = 0
  do ix = 0,NX_GOCAD_HR-2
    do iy = 0,NY_GOCAD_HR-2
      do iz = 0,NZ_GOCAD_HR-2

! suppress elements that are undefined
   if(iflag_point(ix,iy,iz) .and. &
      iflag_point(ix+1,iy,iz) .and. &
      iflag_point(ix+1,iy+1,iz) .and. &
      iflag_point(ix,iy+1,iz) .and. &
      iflag_point(ix,iy,iz+1) .and. &
      iflag_point(ix+1,iy,iz+1) .and. &
      iflag_point(ix+1,iy+1,iz+1) .and. &
      iflag_point(ix,iy+1,iz+1)) then

        ispec = ispec + 1

! use Z > 0 and Z < 0 to define material flag
        zcoord = ORIG_Z_GOCAD_HR + iz*SPACING_Z_GOCAD_HR
        if(zcoord <= 0.) then
          imaterial = 1
        else
          imaterial = 2
        endif

        write(11,200) ispec,imaterial,ipoin_store(ix,iy,iz), &
          ipoin_store(ix+1,iy,iz),ipoin_store(ix+1,iy+1,iz),ipoin_store(ix,iy+1,iz), &
          ipoin_store(ix,iy,iz+1),ipoin_store(ix+1,iy,iz+1), &
          ipoin_store(ix+1,iy+1,iz+1),ipoin_store(ix,iy+1,iz+1)

     endif

      enddo
    enddo
  enddo

 200 format(i6,1x,i2,' hex ',i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)

! output AVS header for data
  write(11,*) '1 1'
  write(11,*) 'Zcoord, meters'

! output data values (P-velocity at points)
  ipoin = 0
  do ix = 0,NX_GOCAD_HR-1
    do iy = 0,NY_GOCAD_HR-1
      do iz = 0,NZ_GOCAD_HR-1
      if(iflag_point(ix,iy,iz)) then
        ipoin = ipoin + 1

! use Vp to color the model
        write(11,*) ipoin,' ',sngl(vp_block_gocad(ix,iy,iz))

! or use Z > 0 and Z < 0 to color the model
!       zcoord = ORIG_Z_GOCAD_HR + iz*SPACING_Z_GOCAD_HR
!       if(zcoord <= 0.) then
!         write(11,*) ipoin,' 0.'
!       else
!         write(11,*) ipoin,' 255.'
!       endif

      endif
      enddo
    enddo
  enddo

  close(11)

  end program read_gocad_block_extract

