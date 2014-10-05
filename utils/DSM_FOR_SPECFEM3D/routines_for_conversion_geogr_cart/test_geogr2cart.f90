  program test_geogr2cart
    use geogr2cart_mod
    implicit none
    double precision x,y,z
    double precision long,lati,depth
    double precision lon_center_chunk,lat_center_chunk, chunk_azi
    double precision ksi,eta

   long=60.
   lati=0.
   depth=0.
   write(*,*) long,lati,depth

   call read_chunk_parameters()
   lon_center_chunk=LON_CENTER
   lat_center_chunk=LAT_CENTER
   chunk_azi=AZI_CHUNK
   call compute_rotation_matrix(rotation_matrix,lon_center_chunk,lat_center_chunk, chunk_azi)
   call  geogr2cart(x,y,z,long,lati,depth,ksi,eta)

    write(*,*) x,y,z

  end program test_geogr2cart
