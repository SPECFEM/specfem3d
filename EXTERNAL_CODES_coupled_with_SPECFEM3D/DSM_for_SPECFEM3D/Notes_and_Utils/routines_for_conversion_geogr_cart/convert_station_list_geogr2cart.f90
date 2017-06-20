  program convert_station_list
    use geogr2cart_mod
    implicit none
    double precision x,y,z
    double precision long,lati,depth
    double precision lon_center_chunk,lat_center_chunk, chunk_azi
    double precision ksi,eta
    character(len=6) cdummy,ST
    integer ista

    1000 format(a6,3x,a2,3x,4f20.5)

    call read_chunk_parameters()
    lon_center_chunk=LON_CENTER
    lat_center_chunk=LAT_CENTER
    chunk_azi=AZI_CHUNK
    call compute_rotation_matrix(rotation_matrix,lon_center_chunk,lat_center_chunk, chunk_azi)

    ista=0
    open(10,file='list_stations_files_to_convert.geogr')
    open(20,file='STATIONS_INPUT_FOR_SEM')
    read(10,'(a)') cdummy
    do
       read(10,*,end=99) lati,long,depth
       x=0.;y=0.;z=0.;ksi=0.;eta=0.
       call  geogr2cart(x,y,z,long,lati,depth,ksi,eta)
       ista=ista+1
       write(ST,'(a1,i5.5)') 'S',ista
       write(20,1000) ST,'TS',x,y,-z,-z
    enddo
99  close(10)
    close(20)
  end program convert_station_list

