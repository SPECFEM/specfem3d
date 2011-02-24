
program testf

implicit none

include 'constants.h'

double precision lon,lat
double precision x,y

integer iway, iproject

iway = ILONGLAT2UTM
iproject = 11

print *, 'input lon, lat'
read(*,*) lon, lat

call utm_geo(lon,lat,x,y,iproject,iway)

print *, 'x = ', x, '  y = ', y

end program testf

