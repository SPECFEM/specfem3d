
program convert_mnt_lacq

!! DK DK convert Lacq MNT from Thomas' format to SPECFEM3D format

implicit none

integer, parameter :: NX_MNT = 499, NY_MNT = 401

double precision, parameter :: DELTAX_MNT = 100.d0, DELTAY_MNT = 100.d0

double precision, parameter :: OFFSETX_MNT = 340400.d0, OFFSETY_MNT = 100000.d0

integer ix,iy

double precision x,y

double precision, dimension(NX_MNT,NY_MNT) :: z,xdiff,ydiff

do ix = 1,NX_MNT
do iy = 1,NY_MNT

read(*,*) x,y,z(ix,iy)

xdiff(ix,iy) = x - ((ix-1)*DELTAX_MNT + OFFSETX_MNT)
ydiff(ix,iy) = y - ((iy-1)*DELTAY_MNT + OFFSETY_MNT)

!!print *,x,sngl(xdiff(ix,iy)),y,sngl(ydiff(ix,iy))

enddo
enddo

!! DK DK write elevation converted to integer to new file
!! DK DK swap ix and iy loops to convert to SPECFEM3D storage convention for MNTs
do iy = 1,NY_MNT
do ix = 1,NX_MNT
  print *,nint(z(ix,iy))
enddo
enddo

print *,'xdiff min max = ',minval(xdiff),maxval(xdiff)
print *,'ydiff min max = ',minval(ydiff),maxval(ydiff)
print *,'z min max = ',minval(z),maxval(z)
print *,'average z = ',sum(z) / dble(NX_MNT*NY_MNT)

end program convert_mnt_lacq

