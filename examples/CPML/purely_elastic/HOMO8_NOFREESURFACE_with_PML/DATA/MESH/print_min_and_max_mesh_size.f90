
 program print_min_and_max_mesh_size

 implicit none

 real, dimension(:), allocatable :: x,y,z

 integer :: ipoin,npoin,idummy

 open(unit=23,file="nodes_coords_file",status="old")
 read(23,*) npoin

 allocate(x(npoin))
 allocate(y(npoin))
 allocate(z(npoin))

 do ipoin = 1,npoin
  read(23,*) idummy,x(ipoin),y(ipoin),z(ipoin)
 enddo

 close(23)

 print *, 'X min max = ',minval(x),maxval(x),' deltaX = ',maxval(x) - minval(x)
 print *, 'Y min max = ',minval(y),maxval(y),' deltaY = ',maxval(y) - minval(y)
 print *, 'Z min max = ',minval(z),maxval(z),' deltaZ = ',maxval(z) - minval(z)

 end program print_min_and_max_mesh_size

