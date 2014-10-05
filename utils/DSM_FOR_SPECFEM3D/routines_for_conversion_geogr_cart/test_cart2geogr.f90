  program test_cart2geogr
    use cart2geogr_mod
    implicit none
    double precision x,y,z
    double precision long,lati,depth
    double precision lon_center_chunk,lat_center_chunk, chunk_azi
    double precision ksi,eta

    x=0.
    y=0.
    z=1101121.76137563



    call read_chunk_parameters()
    lon_center_chunk=LON_CENTER
    lat_center_chunk=LAT_CENTER
    chunk_azi=AZI_CHUNK
    call compute_rotation_matrix(rotation_matrix,lon_center_chunk,lat_center_chunk, chunk_azi)
    call  cart2geogr(x,y,z,long,lati,depth,ksi,eta)


    write(*,*) long,lati,depth

  end program test_cart2geogr
