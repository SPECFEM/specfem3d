program create_interace

  implicit none
  integer i,j,nx,ny
  real hx,hy,x,y,z

  nx=129;ny=29
  hx=25.
  hy=25.
  open(10,file='flat_top.dat')
  do j=1,ny
     do i=1,nx
        z=0.
        write(10,*) z
     end do
  end do
  close(10)
end program create_interace
