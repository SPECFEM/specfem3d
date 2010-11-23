
  program read_gocad_block_extract

  implicit none

  include "../../../constants.h"
  include "../../../constants_gocad.h"

!
! new Voxet Peter July 29, 2002
!

  double precision vp_block_gocad(0:NX_GOCAD_MR-1,0:NY_GOCAD_MR-1,0:NZ_GOCAD_MR-1)
  logical iflag_point(0:NX_GOCAD_MR-1,0:NY_GOCAD_MR-1,0:NZ_GOCAD_MR-1)
  integer ipoin_store(0:NX_GOCAD_MR-1,0:NY_GOCAD_MR-1,0:NZ_GOCAD_MR-1)

! use integer array to store topography values
  integer itopo_bathy_basin(NX_TOPO,NY_TOPO)
  integer iclosestlong,iclosestlat
  double precision elevation,max_error
  double precision lat,long

  integer ix,iy,iz,iz_found,ipoin,ispec,npoin,nspec
  integer icount_undefined

  double precision vpmin,vpmax
  double precision xcoord,ycoord,zcoord
  double precision zsedim_found
  integer irecord,nrecord,i_vp

  print *
  print *,'reading velocity block from Gocad voxet'
  print *

  print *
  print *,'number of points in block NX_GOCAD_MR,NY_GOCAD_MR,NZ_GOCAD_MR = ',NX_GOCAD_MR,NY_GOCAD_MR,NZ_GOCAD_MR
  print *,'total points in block NX_GOCAD_MR*NY_GOCAD_MR*NZ_GOCAD_MR = ',NX_GOCAD_MR*NY_GOCAD_MR*NZ_GOCAD_MR
  print *

! initialize array to undefined values everywhere
  vp_block_gocad(:,:,:) = -1.

! read Vp from extracted text file
  open(unit=27,file='LA_MR_voxet_extracted.txt',status='old')
  read(27,*) nrecord
  do irecord = 1,nrecord
    read(27,*) ix,iy,iz,i_vp
    if(ix<0 .or. ix>NX_GOCAD_MR-1 .or. iy<0 .or. iy>NY_GOCAD_MR-1 .or. iz<0 .or. iz>NZ_GOCAD_MR-1) &
      stop 'wrong array index read in Gocad medium-resolution file'
    vp_block_gocad(ix,iy,iz) = dble(i_vp)
  enddo
  close(27)

!---

  icount_undefined = 0
  vpmin = + 100000000.
  vpmax = - 100000000.
  ipoin = 0

! count total number of points kept
  do ix = 0,NX_GOCAD_MR-1
    do iy = 0,NY_GOCAD_MR-1
      do iz = 0,NZ_GOCAD_MR-1

! exclude points that are undefined
! a negative P velocity has been used to flag these points
        if(vp_block_gocad(ix,iy,iz) < 1.) then
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

! count total number of elements
  ispec = 0
  do ix = 0,NX_GOCAD_MR-2
    do iy = 0,NY_GOCAD_MR-2
      do iz = 0,NZ_GOCAD_MR-2

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
  print *,'which is ',100.*icount_undefined/dble(NX_GOCAD_MR*NY_GOCAD_MR*NZ_GOCAD_MR),' %'
  print *

  print *,'minval maxval Vp in region kept = ',vpmin,vpmax

! create DX file with velocity block
! only save elements that are not fictitious

  print *
  print *,'creating DX file'
  print *,'points = ',npoin
  print *,'elements = ',nspec
  print *

  open(unit=11,file='DX_Pvelocity_block_MR.dx',status='unknown')

  write(11,*) 'object 1 class array type float rank 1 shape 3  items ',npoin,' data follows'

! output global DX points
! loop on all the points
  ipoin = 0
  do ix = 0,NX_GOCAD_MR-1
    do iy = 0,NY_GOCAD_MR-1
      do iz = 0,NZ_GOCAD_MR-1

    if(iflag_point(ix,iy,iz)) then
        ipoin = ipoin + 1

! define Gocad grid, shift of Voxet is taken into account
        xcoord = ORIG_X_GOCAD_MR + ix*SPACING_X_GOCAD_MR
        ycoord = ORIG_Y_GOCAD_MR + iy*SPACING_Y_GOCAD_MR
        zcoord = ORIG_Z_GOCAD_MR + iz*SPACING_Z_GOCAD_MR
        write(11,*) sngl(xcoord),' ',sngl(ycoord),' ',sngl(zcoord)
     endif

      enddo
    enddo
  enddo

  write(11,*) 'object 2 class array type int rank 1 shape 8 items ',nspec,' data follows'

! output global DX elements
! loop on all the elements
  ispec = 0
  do ix = 0,NX_GOCAD_MR-2
    do iy = 0,NY_GOCAD_MR-2
      do iz = 0,NZ_GOCAD_MR-2

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
        zcoord = ORIG_Z_GOCAD_MR + iz*SPACING_Z_GOCAD_MR

! point order in OpenDX is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
        write(11,200) ipoin_store(ix,iy+1,iz)-1,ipoin_store(ix,iy,iz)-1, &
                      ipoin_store(ix,iy+1,iz+1)-1,ipoin_store(ix,iy,iz+1)-1, &
                      ipoin_store(ix+1,iy+1,iz)-1,ipoin_store(ix+1,iy,iz)-1, &
                      ipoin_store(ix+1,iy+1,iz+1)-1,ipoin_store(ix+1,iy,iz+1)-1

     endif

      enddo
    enddo
  enddo

 200 format(i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)

! output DX header for data
  write(11,*) 'attribute "element type" string "cubes"'
  write(11,*) 'attribute "ref" string "positions"'
  write(11,*) 'object 3 class array type float rank 0  items ',npoin,' data follows'

! output data values (P-velocity at points)
  ipoin = 0
  do ix = 0,NX_GOCAD_MR-1
    do iy = 0,NY_GOCAD_MR-1
      do iz = 0,NZ_GOCAD_MR-1
      if(iflag_point(ix,iy,iz)) then
        ipoin = ipoin + 1

! use Vp to color the model
!        write(11,*) sngl(vp_block_gocad(ix,iy,iz))

! or use Z > 0 and Z < 0 to color the model
        zcoord = ORIG_Z_GOCAD_MR + iz*SPACING_Z_GOCAD_MR
        if(zcoord <= 0.) then
          write(11,*) '0'
        else
          write(11,*) '255'
        endif

      endif
      enddo
    enddo
  enddo

  write(11,*) 'attribute "dep" string "positions"'
  write(11,*) 'object "irregular connections  irregular positions" class field'
  write(11,*) 'component "positions" value 1'
  write(11,*) 'component "connections" value 2'
  write(11,*) 'component "data" value 3'
  write(11,*) 'end'

  close(11)

  end program read_gocad_block_extract

