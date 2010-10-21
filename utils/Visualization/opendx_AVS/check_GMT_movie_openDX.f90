
  program check_GMT_movie

! check GMT movie files using OpenDX

  implicit none

  integer, parameter :: NEX_XI = 384+1, NEX_ETA = 384+1

  double precision long(NEX_XI,NEX_ETA)
  double precision lat(NEX_XI,NEX_ETA)
  double precision uz(NEX_XI,NEX_ETA)

  integer ix,iy,iglob1,iglob2,iglob3,iglob4

! read file
  open(unit=13,file='gmt_movie_000001.xyz',status='old')
  do iy=1,NEX_ETA
    do ix=1,NEX_XI
      read(13,*) long(ix,iy),lat(ix,iy),uz(ix,iy)
    enddo
  enddo
  close(13)

! create OpenDX file
  open(unit=3,file='gmt_movie.dx',status='unknown')
  write(3,*) 'object 1 class array type float rank 1 shape 3 items ',NEX_XI*NEX_ETA,' data follows'

  do iy=1,NEX_ETA
    do ix=1,NEX_XI
      write(3,*) sngl(long(ix,iy)),sngl(lat(ix,iy)),' 0'
    enddo
  enddo

  write(3,*) 'object 2 class array type int rank 1 shape 4 items ',(NEX_XI-1)*(NEX_ETA-1),' data follows'

  do iy=1,NEX_ETA-1
    do ix=1,NEX_XI-1
      iglob1 = (iy-1)*NEX_XI + ix
      iglob2 = (iy-1)*NEX_XI + ix+1
      iglob3 = (iy+1-1)*NEX_XI + ix+1
      iglob4 = (iy+1-1)*NEX_XI + ix
      write(3,210) iglob1-1,iglob4-1,iglob2-1,iglob3-1
    enddo
  enddo

  write(3,*) 'attribute "element type" string "quads"'
  write(3,*) 'attribute "ref" string "positions"'
  write(3,*) 'object 3 class array type float rank 0 items ',NEX_XI*NEX_ETA,' data follows'

  do iy=1,NEX_ETA
    do ix=1,NEX_XI
      write(3,*) uz(ix,iy)
    enddo
  enddo

  write(3,*) 'attribute "dep" string "positions"'
  write(3,*) 'object "irregular positions irregular connections" class field'
  write(3,*) 'component "positions" value 1'
  write(3,*) 'component "connections" value 2'
  write(3,*) 'component "data" value 3'
  write(3,*) 'end'

210 format(i6,1x,i6,1x,i6,1x,i6)

  close(3)

  end program check_GMT_movie

