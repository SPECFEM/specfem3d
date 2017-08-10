program create_uniform_sations

   implicit none
   character(len=1) code
   character(len=2) net
   character(len=250) par_file
   real dx,dy,z,xmin,ymin
   integer nx,ny

   real xc,yc,zc
   integer i,j,k
   character(len=5) name_sta

   write(*,*) ' Parameter file name  '
   write(*,*) ' following format :'
   write(*,*)
   write(*,*) ' nx ny       !: (integer) number of stations in x and y direction (m)'
   write(*,*) ' dx dy       !: (real) stations spacing in x and y direction (m)'
   write(*,*) ' xmin, ymin  !: (real) coordinate of the first station (m)'
   write(*,*) ' z           !: (real) depth of all staions (m) '
   write(*,*) ' code        !: (character(len=1)) first letter for station name '
   write(*,*) ' net         !: (character(len=2)) network name '
   write(*,*)
   write(*,*) ' give Parameter file name ? '

   read(*,'(a)') par_file

   open(10,file=trim(par_file))
   read(10,*) nx,ny
   read(10,*) dx,dy
   read(10,*) xmin,ymin
   read(10,*) z
   read(10,'(a)') code
   read(10,'(a)') net
   close(10)

   open(10,file='station_list.txt')

   k = 0
   zc=z
   yc = ymin - dy
   do j = 1, ny
     yc = yc + dy
     xc = xmin - dx
     do i = 1, nx
      k =  k + 1
      xc = xc + dx
       write(name_sta,'(a,i4.4)') code,k
       write(10,'(a5,2x,a2,2x,4f20.5)') name_sta,net,xc,yc,zc,zc
     enddo
   enddo
   close(10)

end program create_uniform_sations
