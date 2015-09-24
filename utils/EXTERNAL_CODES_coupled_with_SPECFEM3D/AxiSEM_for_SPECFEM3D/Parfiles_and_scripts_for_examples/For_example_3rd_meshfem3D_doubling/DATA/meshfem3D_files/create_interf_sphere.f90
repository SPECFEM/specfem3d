program create_interf_files
!
!
! create prem interface on spherical earth to use as input for meshfem3D
! (*.z files)
!

  implicit none
  double precision :: x, y, z, Rt, Rc
  double precision :: R(0:11)
  double precision :: xmin, xmax, ymin, ymax, dx, dy
  integer          :: nx, ny, interf
  character(len=15) :: file_name


!! user parameters ------------------------

  !! sampling (m)
  dx = 5000.d0
  dy = 5000.d0

  !! size  of the domain (m)
  xmin = -500000.d0
  xmax =  500000.d0
  ymin = -500000.d0
  ymax =  500000.d0

!!-----------------------------------------

  !! PREM INTERFACES (radius (m))
   
  R(0)  = 6371000.d0   !! free surface (Earth radius) (m)
  R(1)  = 6356000.d0   !! 1st interface 
  R(2)  = 6346000.d0   !! 2nd interface 
  R(3)  = 6151000.d0   !! 3nd interface 
  R(4)  = 5971000.d0   !! ...
  R(5)  = 5771000.d0
  R(6)  = 5701000.d0
  R(7)  = 5600000.d0
  R(8)  = 3630000.d0
  R(9)  = 3480000.d0
  R(10) = 1221000.d5
  R(11) = 0.d0

  Rt=R(0)
!! creating 2 file for each prem interfaces 
!! *.z is used for meshfem3d *.xyz write the surface file for
!! user purpose : eg :to check or visualize (not used for meshing)

  do interf=0,11
    Rc = R(interf)     

    write(file_name,'(a7,i3.3)') 'Interf_',interf
    open(10,file=trim(file_name)//'.xyz')
    open(11,file=trim(file_name)//'.z')

    y=ymin - dy
    ny = 0
    do while (y <= ymax )
       y = y + dy 
       x = xmin - dx
       ny = ny + 1

       nx = 0
       do while (x <= xmax)
          x = x + dx
          nx=nx+1
          z = sqrt( Rc*Rc - x*x - y*y) - Rt
          write(10,*) x,y,z
          write(11,*) z 
       end do

    end do

    close(10)
    close(11)

  end do

  write(*,*) 
  write(*,*) 'interfaces file created on grid :'
  write(*,*) 'nx ny :', nx, ny
  write(*,*) 'sampling '
  write(*,*) 'dx dy :', dx, dy
  write(*,*) 


end program create_interf_files

