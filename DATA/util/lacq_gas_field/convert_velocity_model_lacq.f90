
program convert_velocity_model_lacq

!! DK DK convert Lacq velocity model from Thomas' format to SPECFEM3D format

implicit none

integer, parameter :: NX_VELOCITY_MODEL_LACQ = 82, NY_VELOCITY_MODEL_LACQ = 64, NZ_VELOCITY_MODEL_LACQ = 150

double precision, parameter :: DELTAX_VELOCITY_MODEL_LACQ = 500.d0, DELTAY_VELOCITY_MODEL_LACQ = 500.d0, DELTAZ_VELOCITY_MODEL_LACQ = 100.d0

double precision, parameter :: OFFSETX_VELOCITY_MODEL_LACQ = 347475.d0, OFFSETY_VELOCITY_MODEL_LACQ = 107725.d0, OFFSETZ_VELOCITY_MODEL_LACQ = 0.d0

integer ix,iy,iz,ival

double precision x,y,z

double precision, dimension(NX_VELOCITY_MODEL_LACQ*NY_VELOCITY_MODEL_LACQ*NZ_VELOCITY_MODEL_LACQ) :: vp,xdiff,ydiff,zdiff


ival = 0

!! DK DK write number of records to new file in SPECFEM3D velocity model format
print *,NX_VELOCITY_MODEL_LACQ*NY_VELOCITY_MODEL_LACQ*NZ_VELOCITY_MODEL_LACQ

!! DK DK invert iz to convert from depth (Thomas) to height (SPECFEM3D)
do iz = NZ_VELOCITY_MODEL_LACQ,1,-1
do iy = 1,NY_VELOCITY_MODEL_LACQ
do ix = 1,NX_VELOCITY_MODEL_LACQ

ival = ival + 1

read(*,*) x,y,z,vp(ival)

xdiff(ival) = x - ((ix-1)*DELTAX_VELOCITY_MODEL_LACQ + OFFSETX_VELOCITY_MODEL_LACQ)
ydiff(ival) = y - ((iy-1)*DELTAY_VELOCITY_MODEL_LACQ + OFFSETY_VELOCITY_MODEL_LACQ)
!!!!!!!!!!zdiff(ival) = z - ((iz-1)*DELTAZ_VELOCITY_MODEL_LACQ + OFFSETZ_VELOCITY_MODEL_LACQ)

!!print *,x,sngl(xdiff(ival)),y,sngl(ydiff(ival)),z,sngl(zdiff(ival))

!! DK DK write final value converted to integer to new file in SPECFEM3D velocity model format
!! DK DK numbering starts at 0 in Gocad, therefore subtract 1 from indices
  print *,ix-1,iy-1,iz-1,nint(vp(ival)*1000.d0)

enddo
enddo
enddo

print *,'xdiff min max = ',minval(xdiff),maxval(xdiff)
print *,'ydiff min max = ',minval(ydiff),maxval(ydiff)
!!!!!!!!!!print *,'zdiff min max = ',minval(zdiff),maxval(zdiff)
print *,'vp min max = ',minval(vp),maxval(vp)

end program convert_velocity_model_lacq

