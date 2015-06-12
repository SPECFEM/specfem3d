program create_station_list_in_geogr
  implicit none
  integer nx,ny
  real dx,dy
  real xmin,ymin
  real r,long,lat
  integer i,j

  r=6371.

  open(10,file='input_sta.txt')
  read(10,*) nx,ny
  read(10,*) dx,dy
  read(10,*) xmin,ymin
  close(10)

  open(10,file='stations_to_convert.txt')
  lat=ymin-dy
  do j=1,ny
     lat = lat + dy
     long=xmin - dx
     do i=1,nx
        long = long + dx
        write(10,*) r, lat, long
     enddo
  enddo

end program create_station_list_in_geogr
