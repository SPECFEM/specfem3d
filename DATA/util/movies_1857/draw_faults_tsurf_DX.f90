
! create DX file to check t-surfs from GoCaD

! compile with " ifc -Vaxlib " to include library for system calls

 program combine_AVS_DX

 implicit none

 integer ntotpoinAVS_DX    ! number of VRTX
 integer ntotspecAVS_DX    ! number of TRGL

 integer ipoin,ispec,idummy
 integer iglob1,iglob2,iglob3

 double precision xval,yval,zval

! filter input GoCad file
 call system("grep VRTX input.ts | sed -e '1,$s/PVRTX//g' | sed -e '1,$s/VRTX//g' > points.dat")
 call system("grep TRGL input.ts | sed -e '1,$s/TRGL//g' > triangles.dat")

! save and read number of points and triangles
 call system("wc -l points.dat > number_points_triangles.dat")
 call system("wc -l triangles.dat >> number_points_triangles.dat")

 open(unit=10,file='number_points_triangles.dat',status='old')
 read(10,*) ntotpoinAVS_DX
 read(10,*) ntotspecAVS_DX
 close(10)

   open(unit=11,file='DX_fullmesh.dx',status='unknown')
!! DK DK also write to ASCII
   open(unit=12,file='ASCII_fullmesh.dat',status='unknown')

! write points
   write(11,*) 'object 1 class array type float rank 1 shape 3 items ',ntotpoinAVS_DX,' data follows'

   write(12,*) ntotpoinAVS_DX
   write(12,*) ntotspecAVS_DX

   open(unit=10,file='points.dat',status='old')

   do ipoin = 1,ntotpoinAVS_DX
       read(10,*) idummy,xval,yval,zval
       write(11,*) xval,yval,zval
       write(12,*) idummy,xval,yval,zval
   enddo

   close(10)

! write elements
   write(11,*) 'object 2 class array type int rank 1 shape 3 items ',ntotspecAVS_DX,' data follows'

   open(unit=10,file='triangles.dat',status='old')

   do ispec = 1,ntotspecAVS_DX
       read(10,*) iglob1,iglob2,iglob3
       write(11,210) iglob1-1,iglob2-1,iglob3-1
       write(12,*) iglob1,iglob2,iglob3
   enddo

   close(10)

 210     format(i6,1x,i6,1x,i6,1x,i6)

   write(11,*) 'attribute "element type" string "triangles"'
   write(11,*) 'attribute "ref" string "positions"'
   write(11,*) 'object 3 class array type float rank 0 items ',ntotspecAVS_DX,' data follows'

! write data
   do ispec = 1,ntotspecAVS_DX
       write(11,*) '0'
   enddo

   write(11,*) 'attribute "dep" string "connections"'
   write(11,*) 'object "irregular positions irregular connections" class field'
   write(11,*) 'component "positions" value 1'
   write(11,*) 'component "connections" value 2'
   write(11,*) 'component "data" value 3'
   write(11,*) 'end'

 close(11)
 close(12)

! suppress temporary files
 call system("rm -f points.dat triangles.dat number_points_triangles.dat")

 end program combine_AVS_DX

