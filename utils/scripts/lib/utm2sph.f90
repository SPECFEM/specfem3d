
program testf

implicit none

include 'constants.h'

double precision lon,lat
double precision x,y

integer iway, iproject

iway = IUTM2LONGLAT
iproject = 11

print *, 'input x, y'
read(*,*) x, y

call utm_geo(lon,lat,x,y,iproject,iway)

print *, 'lon = ', lon, '  lat = ', lat

end program testf

