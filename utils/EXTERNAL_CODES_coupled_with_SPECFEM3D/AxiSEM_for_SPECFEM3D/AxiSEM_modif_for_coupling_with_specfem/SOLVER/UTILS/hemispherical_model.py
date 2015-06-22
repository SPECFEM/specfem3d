#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

# Define layer boundaries (one more than layers)
layers = [1217.5, 1190., 1160., 1100.]
# Define angles of hemispherical boundaries with a linearly interpolated region in between
angles = [[45., 55.], [50., 60.], [55., 65.]]

vp_in = np.ones((len(layers) - 1, len(angles)))

# define perturbations in each layer and hemisphere
vp_in[0,0] = 0.
vp_in[0,1] = 2.

vp_in[1,0] = 0.
vp_in[1,1] = 5.

vp_in[2,0] = 0.
vp_in[2,1] = 2.

# number of points in theta direction (fine sampling usefull in combination with nearest
# neighbour interpolation)
ntheta = 721
dtheta = 180. / (ntheta - 1)

nlayers = len(layers) - 1
# number of radial point per layer. with nearest neighbour interpolation 2 is fine
nlpl = 2
# distance of points from layer boundaries (e.g. to avoid perturbations on both sides of a
# discontinuity)
dr = .01

f = open('model.sph', 'w')

vp = 0.
vs = 0.
rho = 0.

# total number of points. +1 for the additional zero layer at the bottom
npoints = (nlayers * nlpl + 1)  * ntheta
print >> f, npoints

# write model file
for l in np.arange(nlayers):
    for r in np.linspace(layers[l] - dr, layers[l+1] + dr, nlpl):
        for theta in np.linspace(0., 180., ntheta):
            if theta < angles[l][0]:
                vp = vp_in[l,0]
            elif theta > angles[l][1]:
                vp = vp_in[l,1]
            else:
                # linear interpolation in the central region
                vp = vp_in[l,0] \
                   + (vp_in[l,1] - vp_in[l,0]) / (angles[l][1] - angles[l][0]) \
                        * (theta - angles[l][0])

            print >> f, '%7.2f %6.2f %5.2f %5.2f  %5.2f ' % (r, theta, vp, vs, rho)

# additional zero (relative perturbation!) layer at the bottom to make sure the last layer
# does not extent to the next element boundary. Same approach might be usefull for the
# first layer, but in this case it is the ICB anyway
vp = 0.
r = layers[-1] - dr
for theta in np.linspace(0., 180., ntheta):
    print >> f, '%7.2f %6.2f %5.2f %5.2f  %5.2f ' % (r, theta, vp, vs, rho)

f.close
