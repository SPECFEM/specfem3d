program create_topo

  implicit none
  double precision x,y,z,Rt,Rc
  double precision R1,R2,R3
  double precision xmin,xmax,ymin,ymax,dx,dy
  integer nx,ny
  
  Rt=6371000.d0   !! free surface (Earth radius) (m)
  R1=6351000.d0   !! 1st interface 
  R2=6336000.d0   !! 2nd interface 
  R3=62510000.0   !! 3nd interface 
  !6151000.d0   !! interface considered (m)

  !! sampling (m)
  dx=5000.d0   
  dy=5000.d0
  
  !! size  of the domain (m)
  xmin=-500000.d0
  xmax= 500000.d0
  ymin=-500000.d0
  ymax= 500000.d0
  
  Rc = Rt
  open(10,file='topo_shpere.xyz')
  y=ymin - dy
  ny = 0
  do while (y <= ymax )
     y = y + dy 
     x = xmin - dx
     ny = ny + 1
     nx=0
     do while (x <= xmax)
        x = x + dx
        nx=nx+1
        z = sqrt( Rc*Rc - x*x - y*y) - Rt
        write(10,*) x,y,z
     end do
  end do
  close(10)
  write(*,*) 'nx ny :',nx, ny


  Rc = R1
  open(10,file='interf_20km_shpere.xyz')
  y=ymin - dy
  ny = 0
  do while (y <= ymax )
     y = y + dy
     x = xmin - dx
     ny = ny + 1
     nx=0
     do while (x <= xmax)
        x = x + dx
        nx=nx+1
        z = sqrt( Rc*Rc - x*x - y*y) - Rt
        write(10,*) x,y,z
     end do
  end do
  close(10)
  write(*,*) 'nx ny :',nx, ny

  Rc = R2
  open(10,file='interf_Moho_shpere.xyz')
  y=ymin - dy
  ny = 0
  do while (y <= ymax )
     y = y + dy
     x = xmin - dx
     ny = ny + 1
     nx=0
     do while (x <= xmax)
        x = x + dx
        nx=nx+1
        z = sqrt( Rc*Rc - x*x - y*y) - Rt
        write(10,*) x,y,z
     end do
  end do
  close(10)
  write(*,*) 'nx ny :',nx, ny

  Rc = R3
  open(10,file='interf_120km_shpere.xyz')
  y=ymin - dy
  ny = 0
  do while (y <= ymax )
     y = y + dy
     x = xmin - dx
     ny = ny + 1
     nx=0
     do while (x <= xmax)
        x = x + dx
        nx=nx+1
        z = sqrt( Rc*Rc - x*x - y*y) - Rt
        write(10,*) x,y,z
     end do
  end do
  close(10)
  write(*,*) 'nx ny :',nx, ny


end program create_topo

