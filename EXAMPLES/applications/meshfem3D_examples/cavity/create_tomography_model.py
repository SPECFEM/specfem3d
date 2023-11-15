#!/usr/bin/env python
#
# creates a simple tomography file for this example
#
from __future__ import print_function

import sys
import os

# model dimensions
size_x = 144.0  # (m)
size_y = 288.0  # (m)
size_z = 72.0   # (m)

def create_tomo_model(replace=False):
    """
    creates sample tomographic model for SPECFEM3D
    """
    global size_x,size_y,size_z

    # origin points
    ORIG_X = 0.0
    ORIG_Y = 0.0
    ORIG_Z = 0.0

    # end points
    END_X = ORIG_X + size_x
    END_Y = ORIG_Y + size_y
    END_Z = - size_z          # depth in negative z-direction

    # spacing of given tomography points
    SPACING_X = 4.0
    SPACING_Y = 4.0
    SPACING_Z = -2.0

    # number of cell increments
    NX = 36
    NY = 72
    NZ = 36

    # min/max values
    VP_MIN = 2500.0
    VP_MAX = 8500.0
    VS_MIN = 1500.0
    VS_MAX = 7500.0
    RHO_MIN = 1500.0
    RHO_MAX = 1500.0

    # creates directory
    cmd = 'mkdir -p DATA/tomo_files'
    os.system(cmd)

    # opens file
    file = './DATA/tomo_files/tmp.xyz'
    try:
        f = open(file,'w')
    except:
        print("Error opening file ",file)
        sys.tracebacklimit=1
        raise Exception('file does not open: %s' % file)

    # header info
    print("creating header info...")

    f.write("#------------------------\n")
    f.write("# Sample tomographic file\n")
    f.write("#------------------------\n")
    f.write("#orig_x orig_y orig_z end_x end_y end_z\n")
    f.write("%f %f %f %f %f %f\n" % (ORIG_X,ORIG_Y,ORIG_Z,END_X,END_Y,END_Z))
    f.write("#spacing_x spacing_y spacing_z\n")
    f.write("%f %f %f\n" % (SPACING_X,SPACING_Y,SPACING_Z))
    f.write("#nx ny nz\n")
    f.write("%d %d %d\n" % (NX,NY,NZ))
    f.write("#vpmin vpmax vsmin vsmax rhomin rhomax\n")
    f.write("%f %f %f %f %f %f\n" % (VP_MIN,VP_MAX,VS_MIN,VS_MAX,RHO_MIN,RHO_MAX))

    # velocity gradient
    GRADIENT = 0.1
    VS_RATIO = 1.73  # sqrt(3) for Poisson solid

    # adds point location and velocity model values
    print("adding model values...")
    f.write("#model values\n")

    # format: lists first all x, then y, then z
    for k in range(0,NZ):
        for j in range(0,NY):
            for i in range(0,NX):
                x = i * SPACING_X
                y = j * SPACING_Y
                z = k * SPACING_Z
                vp  = VP_MIN + GRADIENT * (-z)
                vs  = VS_MIN + GRADIENT * (-z)/VS_RATIO
                rho = RHO_MIN
                f.write("%f %f %f %f %f %f\n" % (float(x),float(y),float(z),float(vp),float(vs),float(rho)))
    f.close()

    # replaces tomography file
    if replace:
        # renames file
        cmd = 'mv -v DATA/tomo_files/tmp.xyz DATA/tomo_files/tomography_model.xyz'
        os.system(cmd)
        print("created file: ./DATA/tomo_files/tomography_model.xyz")
    else:
        print("see file: ./DATA/tomo_files/tmp.xyz")
    print("")


def usage():
    print("usage:")
    print("    ./create_tomography_model.py <replace>")
    print("  with")
    print("    replace         = flag to force replacing of file [0==check-only/1==replace]")

#
#----------------------------------------------------------------------------
#

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) != 2:
        usage()
        sys.exit(1)
    else:
        if int(sys.argv[1]) == 1:
            replace = True
        else:
            replace = False

    create_tomo_model(replace=replace)
