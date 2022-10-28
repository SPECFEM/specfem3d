#!/usr/bin/env python
import numpy as np

nstat = 91

lat_min = -90
lat_max = 90

lon = np.zeros(nstat)
lat = np.linspace(lat_min, lat_max, nstat)

f = open('STATIONS', 'w')

for i in np.arange(nstat):
    print >> f, 'R%03d   RS      %7.3f   %7.3f   0.0   0.0' % (i, lat[i], lon[i])


f.close()
