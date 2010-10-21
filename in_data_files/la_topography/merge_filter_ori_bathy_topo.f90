
  program merge_filter_ori_bathy_topo

!! DK DK merge original bathy and topo files, truncate and filter result

  implicit none

  include "../../constants.h"

! filter final surface using box filter
  logical, parameter :: FILTER_USING_BOX = .true.
  integer, parameter :: SIZE_FILTER_ONE_SIDE = 20

! truncate high mountains in Sierra Nevada
  logical, parameter :: TRUNCATE_TOPO = .true.

! size of bathymetry grid
  integer, parameter :: NLINES_BATHY_FILE = 29462
  integer, parameter :: NX_BATHY = 241,NY_BATHY = 241
  double precision, parameter :: DEGREES_PER_CELL_BATHY = 5.d0/300.d0

  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_read,ibathy
  integer, dimension(NX_TOPO,NY_TOPO) :: itopo,itopo_filtered

  integer ix,iy,ixconv,iyconv,i,ic,ixcount,iycount
  integer ixbox,iybox,ixval,iyval

  double precision rlon,rlat,a,b,sumval

! read bathy
  ibathy_read(:,:) = 0
  ibathy(:,:) = 0
  print *,'reading ori bathy file'
  open(unit=13,file='bathy_LA_2745618chr_ori.dat',status='old')
  do i=1,NLINES_BATHY_FILE
    read(13,*) a,b,ic
    ix = nint((b-ORIG_LONG)/DEGREES_PER_CELL_BATHY)
    iy = nint((a-ORIG_LAT)/DEGREES_PER_CELL_BATHY)
    if(ix < 1) ix = 1
    if(ix > NX_BATHY) ix = NX_BATHY
    if(iy < 1) iy = 1
    if(iy > NY_BATHY) iy = NY_BATHY
    if(ic <= 0) stop 'incorrect bathy point'
    ibathy_read(ix,iy) = - ic
  enddo
  close(13)

  ibathy = ibathy_read

! remove zeros (spikes) from raw file
  do iy=1,NY_BATHY
    do ix=1,NX_BATHY-1
      if(ibathy_read(ix,iy) >= 0) ibathy(ix,iy) = ibathy_read(ix+1,iy)
    enddo
  enddo

  do iy=1,NY_BATHY-1
    do ix=1,NX_BATHY
      if(ibathy_read(ix,iy) >= 0) ibathy(ix,iy) = ibathy_read(ix,iy+1)
    enddo
  enddo

! remove first column, which has artefacts
  do iy=1,NY_BATHY
    ibathy(1,iy) = ibathy(5,iy)
    ibathy(2,iy) = ibathy(5,iy)
    ibathy(3,iy) = ibathy(5,iy)
    ibathy(4,iy) = ibathy(5,iy)
  enddo

! remove first row, which has artefacts
  do ix=1,NX_BATHY
    ibathy(ix,1) = ibathy(ix,5)
    ibathy(ix,2) = ibathy(ix,5)
    ibathy(ix,3) = ibathy(ix,5)
    ibathy(ix,4) = ibathy(ix,5)
  enddo


! read topo
  print *,'reading ori topo file'
  open(unit=13,file='topo_18sec_121_114_32_37_ori.dat',status='old')
  do iy = NY_TOPO,1,-1
    do ix = 1,NX_TOPO
      read(13,*) itopo(ix,iy)
    enddo
  enddo
  close(13)

! get interpolated bathy where topo is zero
  print *,'interpolating bathy onto topo grid'
  do iy = 1,NY_TOPO
  do ix = 1,NX_TOPO
    if(itopo(ix,iy) <= 0) then
      rlon = ORIG_LONG + (ix-1)*DEGREES_PER_CELL
      rlat = ORIG_LAT + (iy-1)*DEGREES_PER_CELL
      ixconv = nint((rlon - ORIG_LONG)/DEGREES_PER_CELL_BATHY + 1)
      iyconv = nint((rlat - ORIG_LAT)/DEGREES_PER_CELL_BATHY + 1)
      if(ixconv < 1) ixconv = 1
      if(iyconv < 1) iyconv = 1
      if(ixconv > NX_BATHY) ixconv = NX_BATHY
      if(iyconv > NY_BATHY) iyconv = NY_BATHY
      itopo(ix,iy) = ibathy(ixconv,iyconv)
    endif
  enddo
  enddo

! truncate topo in Sierra Nevada to avoid artefacts
  if(TRUNCATE_TOPO) then
    print *,'truncating topo above ',MAX_TOPO
    do iy = 1,NY_TOPO
      do ix = 1,NX_TOPO
        if(itopo(ix,iy) > MAX_TOPO) itopo(ix,iy) = MAX_TOPO
      enddo
    enddo
  endif

! filter final surface using box filter
  if(FILTER_USING_BOX) then
    print *,'filtering final surface using box filter'
    do iy = 1,NY_TOPO
      print *,'doing iy = ',iy,' out of ',NY_TOPO
      do ix = 1,NX_TOPO
        sumval = 0.d0
        do iybox = iy-SIZE_FILTER_ONE_SIDE,iy+SIZE_FILTER_ONE_SIDE
          do ixbox = ix-SIZE_FILTER_ONE_SIDE,ix+SIZE_FILTER_ONE_SIDE
            ixval = ixbox
            iyval = iybox
            if(ixval < 1) ixval = 1
            if(iyval < 1) iyval = 1
            if(ixval > NX_TOPO) ixval = NX_TOPO
            if(iyval > NY_TOPO) iyval = NY_TOPO
            sumval = sumval + dble(itopo(ixval,iyval))
          enddo
        enddo
        itopo_filtered(ix,iy) = nint(sumval/dble((2*SIZE_FILTER_ONE_SIDE+1)**2))
      enddo
    enddo
    itopo(:,:) = itopo_filtered(:,:)
  endif

! output final merged file
  print *,'saving final merged file'
  open(unit=13,file='topo_bathy_final.dat',status='unknown')
  do iy = 1,NY_TOPO
    do ix = 1,NX_TOPO
      write(13,*) itopo(ix,iy)
    enddo
  enddo
  close(13)

! output final merged file
  print *,'saving subsampled merged file for OpenDX'
  open(unit=13,file='topo_bathy_subsampled_OpenDX.dat',status='unknown')
  iycount = 0
  do iy = 1,NY_TOPO,SUBSAMP_FACTOR_OPENDX
    iycount = iycount + 1
    ixcount = 0
    do ix = 1,NX_TOPO,SUBSAMP_FACTOR_OPENDX
      ixcount = ixcount + 1
      write(13,*) itopo(ix,iy)
    enddo
  enddo
  close(13)
  print *,'size of subsampled file for OpenDX is ',ixcount,iycount

  end program merge_filter_ori_bathy_topo

