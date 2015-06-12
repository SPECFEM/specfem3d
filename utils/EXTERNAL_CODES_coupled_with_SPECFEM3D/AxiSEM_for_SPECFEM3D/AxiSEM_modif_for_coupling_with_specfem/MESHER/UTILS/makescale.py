#!/usr/bin/env python
import numpy as np


nproc = 1104
nradialslices = 12

nrep = nproc / nradialslices / 2

f = open('hot_scale_%04d.xml' % (nrep,), 'w')

print >> f, '<ColorMap name="veryhot" space="Lab" indexedLookup="false">'
for i in np.arange(nrep):
    x = i + 0.
    print >> f, '    <Point x="%5.3f" o="0" r="0" g="0" b="0"/>' % (x,)
    x = i + 1. / 6.
    print >> f, '    <Point x="%5.3f" o="0.4" r="0.901961" g="0" b="0"/>' % (x,)
    x = i + 2. / 6.
    print >> f, '    <Point x="%5.3f" o="0.8" r="0.901961" g="0.901961" b="0"/>' % (x,)
    x = i + 3. / 6.
    print >> f, '    <Point x="%5.3f" o="1" r="1" g="1" b="1"/>' % (x,)
    x = i + 4. / 6.
    print >> f, '    <Point x="%5.3f" o="0.8" r="0.901961" g="0.901961" b="0"/>' % (x,)
    x = i + 5. / 6.
    print >> f, '    <Point x="%5.3f" o="0.4" r="0.901961" g="0" b="0"/>' % (x,)

print >> f, '<NaN r="0" g="0.498039" b="1"/>'

print >> f, '</ColorMap>'

f.close()
