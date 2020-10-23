#!/usr/bin/env python
#
# takes as input a CMTSOLUTION file and plots a beachball representation
#
# required modules:
#   - obspy : install for example via pip install -U obspy
#
from __future__ import print_function

import sys
import os
import datetime

from math import sin,asin,cos,acos,atan2,sqrt,log10,pi

import numpy as np
from numpy.linalg import eigh

try:
    from obspy.imaging.beachball import beachball
except:
    print("importing module obspy.imaging.beachball failed, please make sure to install obspy, for example via: pip install -U obspy")
    sys.exit(1)
#from obspy import read_events


#--------------------------------------------------------------------------------------------------

def read_CMT(cmt_file):
    """
    reads in CMT solution
    """
    # default CMTSOLUTION format:
    #
    # PDE  1999 01 01 00 00 00.00  5000 2000 -1000 2.2 2.2 test
    # event name:       test
    # time shift:       0.0000
    # half duration:    1.0
    # latorUTM:       5000.0
    # longorUTM:      2000.0
    # depth:          1.0
    # Mrr:       1.000000e+18
    # Mtt:       1.000000e+18
    # Mpp:       1.000000e+18
    # Mrt:      -2.500000e+19
    # Mrp:       4.000000e+17
    # Mtp:      -2.500000e+18

    print("CMT file name: ",cmt_file)

    # checks if file exists
    if not os.path.isfile(cmt_file):
        print("file %s does not exist, please check..." % cmt_file)
        sys.exit(1)

    # opens file
    try:
        f = open(cmt_file,'r')
    except:
        print("error opening file ",cmt_file)
        sys.exit(1)

    raw = ""

    # header (first line)
    header = f.readline()
    raw += header

    # counts file lines
    i = 1
    for line in f:
        raw += line
        i += 1

    #debug
    print("\nraw content:\n",raw)

    # checks file length
    if i < 13:
        print("number of filelines %i too small for a CMTSOLUTION file, exiting..." % i)
        sys.exit(1)
    if i%13 != 0:
        print("number of filelines %i doesn't match for a CMTSOLUTION file, exiting..." % i)
        sys.exit(1)

    number_of_cmts = i / 13
    print("number of CMTs contained in file: ", number_of_cmts)
    print("")

    # header information (taken from first cmt)
    # header format: CMT 2011  8 23 17 51  7.30  37.8400  -77.9600  12.0 5.8 5.8 C201108231751A
    pre = header[0:3]
    info = header[4:len(header)].split() # splits on whitespace

    # checks number of items
    if len(info) != 12:
        print("error CMT header line:")
        print("  pretext info: ",pre)
        print("  header info :",info)
        print("  number of info items: ",len(info))
        sys.exit(1)

    # reads in header, date and time
    year = int(info[0])
    month = int(info[1])
    day = int(info[2])
    hour = int(info[3])
    minutes = int(info[4])
    seconds = float(info[5])
    # location
    lat = float(info[6])
    lon = float(info[7])
    depth = float(info[8])
    # magnitudes
    Mb = float(info[9])
    Mw = float(info[10])
    # identifier
    eventname = str(info[11])

    # calcuates milliseconds and microseconds
    sec = int(seconds)
    millisec = int(seconds * 10**3 - sec * 10**3)
    microsec = int(seconds * 10**6 - sec * 10**6)

    d = datetime.datetime(year, month, day, hour, minutes, sec, microsec)
    time = '{:%Y-%m-%dT%H:%M:%S.%f}'.format(d)

    print("header information:")
    print("  date      : ",time)
    print("  position  : lat/lon/depth = ",lat,"/",lon,"/",depth)
    print("  magnitudes: Mb/Mw = ",Mb,Mw)
    print("  event name: ",eventname)
    print("")

    ## gets moment tensors
    # repositions to beginning
    f.seek(0,0)

    moment_tensor = []
    for i in range(number_of_cmts):
        # PDE  1999 01 01 00 00 00.00  5000 2000 -1000 2.2 2.2 test
        # event name:       test
        # time shift:       0.0000
        # half duration:    1.0
        # latorUTM:       5000.0
        # longorUTM:      2000.0
        # depth:          1.0

        # skip first line
        line = f.readline()
        # skip name
        line = f.readline()
        # skip shift
        line = f.readline()
        # skip duration
        line = f.readline()
        # skip lat
        line = f.readline()
        # skip lon
        line = f.readline()
        # skip depth
        line = f.readline()

        # Mrr:       1.000000e+18
        # Mtt:       1.000000e+18
        # Mpp:       1.000000e+18
        # Mrt:      -2.500000e+19
        # Mrp:       4.000000e+17
        # Mtp:      -2.500000e+18
        mt = [0.0,0.0,0.0,0.0,0.0,0.0]
        for j in range(0,6):
            line = f.readline()
            arr = line.split()
            #print(arr)
            # Harvard convention
            # format: [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]
            mt[j] = float(arr[1])

        moment_tensor.append(mt)

    # closes file
    f.close()

    return moment_tensor


#--------------------------------------------------------------------------------------------------


def get_scalar_moment(mt):
    """
    computes scalar moment M0 in dyne-cm (given moment tensor elements in dyne-cm)
    """
    #
    # the euclidean matrix norm is invariant under rotation.
    # thus, input can be:
    #   Mxx,Myy,Mzz,Mxy,Mxz,Myz
    # or
    #   Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
    #
    # euclidean (or Frobenius) norm of a matrix: M0**2 = sum( Mij**2)
    Mxx = mt[0]
    Myy = mt[1]
    Mzz = mt[2]
    Mxy = mt[3]
    Mxz = mt[4]
    Myz = mt[5]
    scalar_moment = Mxx**2 + Myy**2 + Mzz**2 + 2.0 * Mxy**2 + 2.0 * Mxz**2 + 2.0 * Myz**2

    # adds 1/2 to be coherent with double couple or point sources
    scalar_moment = sqrt(0.5 * scalar_moment)

    # note: scale factor for the moment tensor
    #
    # CMTSOLUTION moment-tensor elements are given in dyne-cm
    # (from Global CMT project, Dziewonski 1981, Ekstrom et al. 2012)
    #
    # 1 dyne is 1 gram * 1 cm / (1 second)^2
    # 1 Newton is 1 kg * 1 m / (1 second)^2
    # thus 1 Newton = 100,000 dynes
    # therefore 1 dyne.cm = 1e-7 Newton.m
    #       and 1 Newton.m = 1e7 dyne.cm

    # return value (in dyne-cm)
    M0 = scalar_moment

    return M0


#--------------------------------------------------------------------------------------------------

def get_moment_magnitude(mt):
    """
    computes the moment magnitude Mw (following USGS convention)
    """
    # scalar moment (in dyne-cm)
    M0 = get_scalar_moment(mt)

    # moment magnitude by Hanks & Kanamori, 1979
    # Mw = 2/3 log( M0 ) - 10.7       (dyne-cm)
    #
    # alternative forms:
    # Mw = 2/3 ( log( M0 ) - 16.1 )   (N-m) "moment magnitude" by Hanks & Kanamori(1979) or "energy magnitude" by Kanamori (1977)
    #
    # Aki & Richards ("Quantitative Seismology",2002):
    # Mw = 2/3 ( log( M0 ) - 9.1 )    (N-m)
    #
    # conversion: dyne-cm = 10**-7 N-m
    #
    # we follow here the USGS magnitude policy:
    # "All USGS statements of moment magnitude should use M = (log M0)/1.5-10.7
    #  for converting from scalar moment M0 to moment magnitude. (..)"
    # see: http://earthquake.usgs.gov/aboutus/docs/020204mag_policy.php

    # this is to ensure M0>0.0 inorder to avoid arithmetic error with log-function
    if M0 > 0.0:
        Mw = 2.0/3.0 * log10(M0) - 10.7
    else:
        # dummy value
        Mw = -99

    return Mw


#--------------------------------------------------------------------------------------------------

def sph2cart(azimuth,elevation,r):
    """
    converts spherical coordinates to cartesian.
    """
    # definition identical to matlab function: https://www.mathworks.com/help/matlab/ref/sph2cart.html
    #
    # radius r, azimuth [0,2pi] and elevation [-pi/2,pi/2] (in rad)

    # note: inputs can be numpy vectors, thus using numpy-functions
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)

    return x,y,z

#--------------------------------------------------------------------------------------------------


def cart2sph(x,y,z):
    """
    converts cartesian coordinates to spherical
    """
    # definition identical to matlab function: https://www.mathworks.com/help/matlab/ref/cart2sph.html
    #
    # returns azimuth,elevation,r
    azimuth = np.arctan2(y,x)
    elevation = np.arctan2(z,np.sqrt(x**2 + y**2))
    r = np.sqrt(x**2 + y**2 + z**2)

    # elevation range [-pi/2,pi/2]
    # azimuth range [0,360]

    return azimuth,elevation,r

#--------------------------------------------------------------------------------------------------


def convert_SDR_to_MT(strike,dip,rake,M0=1.0):
    """
    converts strike,dip,rake to moment tensor
    """
    # convert strike-dip-rake to Harvard-convention moment tensor mt
    # f = strike, d = dip, l = rake
    #
    # 1 Mrr =  Mzz =  Mo sin(2d) sin(l)
    # 2 Mtt =  Mxx = -Mo(sin(d) cos(l) sin(2f) +     sin(2d) sin(l) (sin(f))^2 )
    # 3 Mpp =  Myy =  Mo(sin(d) cos(l) sin(2f) -     sin(2d) sin(l) (cos(f))^2 )
    # 4 Mrt =  Mxz = -Mo(cos(d) cos(l) cos(f)  +     cos(2d) sin(l) sin(f) )
    # 5 Mrp = -Myz =  Mo(cos(d) cos(l) sin(f)  -     cos(2d) sin(l) cos(f) )
    # 6 Mtp = -Mxy = -Mo(sin(d) cos(l) cos(2f) + 0.5 sin(2d) sin(l) sin(2f) )
    #
    # see, e.g.: https://github.com/g2e/seizmo/blob/master/cmt/sdr2mt.m

    # scales with scalar moment M0 (dyne-cm)
    #
    Mrr = M0 * sin(2.0*dip) * sin(rake)
    Mtt = -M0 * ( sin(dip) * cos(rake) * sin(2.0*strike) + sin(2.0*dip) * sin(rake) * sin(strike)**2 )
    Mpp = M0 * ( sin(dip) * cos(rake) * sin(2.0*strike) - sin(2.0*dip) * sin(rake) * cos(strike)**2 )
    Mrt = -M0 * ( cos(dip) * cos(rake) * cos(strike)   + cos(2.0*dip) * sin(rake) * sin(strike) )
    Mrp = M0 * ( cos(dip) * cos(rake) * sin(strike)   - cos(2.0*dip) * sin(rake) * cos(strike) )
    Mtp = -M0 * ( sin(dip) * cos(rake) * cos(2.0*strike) + 0.5 * sin(2.0*dip) * sin(rake) * sin(2.0*strike) )

    # harvard convention
    # format: [Mrr Mtt Mpp Mrt Mrp Mtp]
    mt = [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]

    return mt


#--------------------------------------------------------------------------------------------------



def decompose_moment_tensor(mt,verbose=False):
    """
    decomposes moment-tensor into maximum double-couple + clvd
    """
    # returns the maximum double couple with a compensated linear-vector dipole (clvd) component for the
    # moment tensor(s) in MT.
    #
    # Note that this decomposition is not unique.
    # one could also decompose into a maximum clvd with some double-couple component.
    #
    # here, we follow the decomposition used by the GlobalCMT & USGS groups for defining
    # the fault planes associated with their moment tensors.
    #
    # follows original implementation from seizmo: [dblcpl,clvd,vec] = mt_decomp(mt,'maxdc');
    # see: https://github.com/g2e/seizmo/tree/master/cmt

    # mt = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp ]
    # notation: M11,M22,M33,M12,M13,M23 and Mij = Mji
    M11 = mt[0]
    M22 = mt[1]
    M33 = mt[2]
    M12 = mt[3]
    M13 = mt[4]
    M23 = mt[5]

    # diagonalizes moment tensor using eigenvalues
    # symmetric 3x3 moment-tensor
    mat = np.array([ [M11,M12,M13], [M12,M22,M23], [M13,M23,M33] ])

    # note: eigh-function will list largest negative eigenvalue first, then intermediate,then largest positive
    eigval,eigvec = eigh(mat)

    # diagonal form with eigenvalues
    M1 = eigval[0]
    M2 = eigval[1]
    M3 = eigval[2]
    mt_diag = np.array([M1,M2,M3,0.0,0.0,0.0])

    # isotropic part
    #
    #                tr(M)  0     0
    # M_iso = 1/3 [    0   tr(M)  0     ]
    #                  0    0    tr(M)
    #
    # note: trace is invariant under rotation, thus equivalent to:
    # > trace = mt_diag[0] + mt_diag[1] + mt_diag[2]
    trace =  mt[0] + mt[1] + mt[2]

    # mt = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp ] (or [Mzz, Mxx, Myy, Mxz, -Myz, -Mxy])
    mt_iso = 1.0/3.0 * np.array([trace, trace, trace,0.0,0.0,0.0])

    # deviatoric part
    #
    #   M = M_iso + M_dev
    #
    mt_dev = np.array([mt[0],mt[1],mt[2],mt[3],mt[4],mt[5]]) - mt_iso

    # double-couple part
    #
    # decomposition:          isotropic              + double couple                   +  CLVD
    #
    #
    #       M1 0  0            tr(M)  0    0                              0   0  0               -M3  0    0
    # M = [ 0  M2 0  ] = 1/3 [  0   tr(M)  0    ]    + (1 - 2 * eps) * [  0  -M3 0  ]  + eps * [  0  -M3   0  ]
    #       0  0  M3            0     0   tr(M)                           0   0  M3               0   0  2 M3
    #
    # with eps strength of CLVD component (between 0 and 0.5)

    # eigenvalues of deviatoric: M1 - 1/3 tr(M), M2 - 1/3 tr(M), M3 - 1/3 tr(M)
    # also compare with:
    # > mat = np.array([ [mt_dev[0],mt_dev[3],mt_dev[4]], [mt_dev[3],mt_dev[1],mt_dev[5]], [mt_dev[4],mt_dev[5],mt_dev[2]] ])
    # > eigval,eigvec = eigh(mat)
    # > print(mt_diag - mt_iso)
    # > print(eigval)
    M1_dev = M1 - mt_iso[0]
    M2_dev = M2 - mt_iso[1]
    M3_dev = M3 - mt_iso[2]

    if verbose:
        print("  trace = ",trace)
        print("  eigenvalues :",M1,M2,M3)
        print("  eigenvalues deviatoric: ",M1_dev,M2_dev,M3_dev)
        print("")

    # dziewonski et al 1981:
    #   used by globalcmt, USGS catalogs (T, P, B axes shared but mixed)
    #   best double couple is 1/2 the difference of largest pos & neg eigenvalues
    #   clvd is 2 min(abs(egn))/max(abs(egn))

    # max/min eigenvalues
    max_v = max(M1_dev,M2_dev,M3_dev)
    min_v = min(M1_dev,M2_dev,M3_dev)

    # index of max/min eigenvalue
    max_index = mt_diag.argmax()
    min_index = mt_diag.argmin()

    eps = 0.5 * (max_v - min_v)

    # double-couple moment-tensor (in diagonal form)
    mt_DC = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
    mt_DC[max_index] = eps
    mt_DC[min_index] = -eps

    # compensated linear vector dipole (in diagonal form)
    mt_CLVD = mt_diag - mt_iso - mt_DC

    # note: mt_dev is in original axis-system,
    #       mt_DC and mt_CLVD are in diagonal form, which relates to eigenvector-system (eigvec)

    # rotate diagonal moment-tensor back into harvard convention using eigenvectors
    mt_DC_rot = rotate_moment_tensor(mt_DC,eigvec)
    mt_CLVD_rot = rotate_moment_tensor(mt_CLVD,eigvec)

    # note: due to round-off errors, the rotated double-couple matrix might have a non-zero trace
    # imposes a zero-trace constraint for double-couple solution
    trace = mt_DC_rot[0] + mt_DC_rot[1] + mt_DC_rot[2]
    mt_DC_rot[0] -= 1.0/3.0 * trace
    mt_DC_rot[1] -= 1.0/3.0 * trace
    mt_DC_rot[2] -= 1.0/3.0 * trace

    #print("trace DC original : ",trace)
    #print("trace DC corrected: ",mt_DC_rot[0]+mt_DC_rot[1]+mt_DC_rot[2])

    if verbose:
        print("moment tensor decomposition: (Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)")
        print("  MT isotropic    :\n",mt_iso)
        print("  MT deviatoric   :\n",mt_dev)
        #print("")
        #print("  diagonal form:")
        #print("  MT double-couple:\n",mt_DC)
        #print("  MT CLVD         :\n",mt_CLVD)
        #print("")
        #print("  Harvard convention:")
        print("  MT DC           :\n",mt_DC_rot)
        print("  MT CLVD         :\n",mt_CLVD_rot)
        print("")

    return mt_iso,mt_dev,mt_DC_rot,mt_CLVD_rot

#--------------------------------------------------------------------------------------------------

def rotate_moment_tensor(mt,R):
    """
    rotates moment-tensor mt using rotation matrix R
    """
    # converts from eigenvectors & eigenvalues representation to harvard-convention
    # mt = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp ] (given in diagonal-form)
    M11 = mt[0]
    M22 = mt[1]
    M33 = mt[2]
    M12 = mt[3]
    M13 = mt[4]
    M23 = mt[5]

    # symmetric 3x3 moment-tensor
    mat = np.array([ [M11,M12,M13], [M12,M22,M23], [M13,M23,M33] ])

    # using rotation matrix R (is eigenvector matrix here):  R * M * transp(R)
    res = np.matmul(mat,np.transpose(R))
    mat = np.matmul(R,res)

    # note: rotation back to original coordinate-frame will introduce numerical round-off errors
    #print("  rotated mat:\n",R,mat)

    # mt = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp ]
    mt_rot = np.array([mat[0,0],mat[1,1],mat[2,2],mat[0,1],mat[0,2],mat[1,2]])

    return mt_rot

#--------------------------------------------------------------------------------------------------

def principal_axis(mt,verbose=False):
    """
    returns principal axis of moment-tensor mt
    """
    # principle axis
    # t, n, p or T, P, B (the tension, pressure & null axes)
    # characterized by eigenvalue, plunge and azimuth, where plunge & azimuth are in degrees:
    # plunge : positive downward from the horizontal
    # azimuth: positive clockwise from North
    #
    # follows: https://github.com/g2e/seizmo/blob/master/cmt/mt2tpb.m
    #         [t,p,b] = mt2tpb(mt_undiag(dblcpl,vec));

    # [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp ]
    M11 = mt[0]
    M22 = mt[1]
    M33 = mt[2]
    M12 = mt[3]
    M13 = mt[4]
    M23 = mt[5]

    # symmetric 3x3 moment-tensor
    mat = np.array([ [M11,M12,M13], [M12,M22,M23], [M13,M23,M33] ])

    # returns eigenvalue,eigenvectors
    # note: eigenvalues are normalized
    d,vec = eigh(mat)

    # convert eigenvectors in Harvard format to plunge and azimuth
    x = vec[1]
    y = -vec[2]
    z = vec[0]

    az,ele,r = cart2sph(x,y,z)

    # angles in degree
    az = az * 180.0/pi      # azimuth in degrees
    pl = ele * 180.0/pi     # elevation -> plunge

    #print("x,y,z:",x,y,z)
    #print("az,pl:",az,pl)

    # round-off errors might lead to small variations for plunge around zero. (for z-coordinate close to zero)
    # this will affect the azimuth (+/- 180)
    #
    # convention plunge always positive, azimuth in [0,360]
    # rotate azimuth 180 if plunge is negative
    for i in range(0, 3):
        if pl[i] < 0.0:
            pl[i] = -pl[i]
            az[i] += 180.0
        if az[i] < 0.0:
            az[i] += 360.0
        if az[i] > 360.0:
            az[i] -= 360.0

    # note: numpy eigh-function lists largest negative eigenvalue first, then intermediate, then largest positive
    #
    # vector format: value, azimuth, plunge
    t = [d[2], az[2], pl[2]]  # tension axis     <-> largest positive eigenvalue
    n = [d[1], az[1], pl[1]]  # null axis        <-> "intermediate" eigenvalue
    p = [d[0], az[0], pl[0]]  # compression axis <-> largest negative eigenvalue

    ## compare with obspy version
    ## see function mt2axes: https://github.com/obspy/obspy/blob/master/obspy/imaging/beachball.py
    ## should give identical results...
    #print("t n p:\n",t,n,p)
    #az = np.arctan2(vec[2], -vec[1])
    #pl = np.arcsin(-vec[0])
    #for i in range(0, 3):
    #    if pl[i] <= 0.0:
    #        pl[i] = -pl[i]
    #        az[i] += pi
    #    if az[i] < 0.0:
    #        az[i] += 2.0 * pi
    #    if az[i] > 2.0 * pi:
    #        az[i] -= 2.0 * pi
    #pl *= 180.0/pi
    #az *= 180.0/pi
    #t = [d[2], az[2], pl[2]]
    #n = [d[1], az[1], pl[1]]
    #p = [d[0], az[0], pl[0]]
    #print("t n p obspy:\n",t,n,p)

    if verbose:
        print("principal axes: (eigenvalue, azimuth, plunge)")
        print("  tension     t = ",t)
        print("  compression p = ",p)
        print("  null        n = ",n)
        print("")

    # checks if valid
    if t[0] == n[0] and n[0] == p[0] and t[0] == 0.0:
        print("WARNING: all eigenvalues are zero! there is no double-couple mechanism involved in this moment tensor.\n")
        return t,n,p

    if abs(pl[0]) < 1.e-5:
        print("WARNING: plunge of compression axis is close to zero! this will affect the azimuth by +/- 180 degrees\n")
    if abs(pl[2]) < 1.e-5:
        print("WARNING: plunge of tension axis is close to zero! this will affect the azimuth by +/- 180 degrees\n")

    return t,n,p


#--------------------------------------------------------------------------------------------------

def fault_plane(t,n,p,verbose=False):
    """
    determines normal and slip vector for fault plane, given the principal axis t/n/p (tension,null,compression)
    """
    # tension axis
    # converts vector in [Value/Azimuth/Plunge] to [North/East/Up]
    v = t[0]      # eigenvalue
    az = t[1]     # azimuth
    pl = t[2]     # plunge

    # tension axis in North/East/Up
    r = 1.0                       # takes unit vector
    az = - az * pi/180.0 + pi/2.0 # azimuth
    ele = - pl * pi/180.0         # elevation/plunge

    x,y,z = sph2cart(az,ele,r)  # x == east / y == north / z == up
    t_neu = np.array([y,x,z])   # north/east/up

    # pressure axis
    v = p[0]
    az = p[1]
    pl = p[2]

    # pressure axis in North/East/Up
    r = 1.0                       # takes unit vector
    az = - az * pi/180.0 + pi/2.0 # azimuth
    ele = - pl * pi/180.0         # elevation/plunge

    x,y,z = sph2cart(az,ele,r)
    p_neu = np.array([y,x,z])   # north/east/up

    #print("t:",t_neu)
    #print("p:",p_neu)

    # gets normal & slip (for fault plane 1)
    normal = t_neu + p_neu        # for plane 2: t_neu - p_neu
    slip   = t_neu - p_neu        # for plane 2: t_neu + p_neu

    # precision fix
    #
    # Converting between t & p axes of vertical faults (dip=90deg)
    # fault normal can give a slightly downwards normal
    #if normal[2] < 0.0 and normal[2] > -1.e-31: normal[2] = 0.0

    # if normal is pointing downwards, flip the sign of both the normal and
    # slip vectors to give the appropriate normal and slip vectors (those corresponding to the hanging wall)
    if normal[2] < 0.0:
        normal = - normal
        slip = - slip

    if verbose:
        print("fault plane:")
        print("  normal vector = ",normal)
        print("  slip vector   = ",slip)
        print("")

    # normalize slip vector
    n = sqrt(slip[0]**2 + slip[1]**2 + slip[2]**2)
    if n > 0.0:
        slip = slip / n

    # normalize normal vector
    n = sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)
    if n > 0.0:
        normal = normal / n

    # now: normal is the normal vector to the fault and must point to the hanging wall (thus z-component is always >= 0)
    if normal[2] < 0.0:
        print("error: invalid normal, must point up: ",normal)
        sys.exit(1)

    return normal,slip


#--------------------------------------------------------------------------------------------------


def convert_principal_axis_to_SDR(t,n,p,verbose=False):
    """
    computes strike,dip,rake from principal axis t (tension), n (null), p (compression/pressure) axis
    """
    # strike: clockwise from North
    # dip   : positive downward from the horizontal
    # rake  : counter-clockwise in the fault plane from the strike direction
    #
    # Note that the strike must be such that when you look along the direction of the strike the fault dips to your right.
    #[strike,dip,rake] = tpb2sdr(t,p,b);
    #
    # follows: https://github.com/g2e/seizmo/blob/master/cmt/tpb2sdr.m

    # checks if double-couple valid
    v1 = t[0]     # eigenvalue
    v2 = n[0]
    v3 = p[0]
    if v1 == v2 and v2 == v3 and v1 == 0.0:
        strike = 0.0
        dip = 0.0
        rake = 0.0
        return strike,dip,rake

    # gets normal and slip vectors for fault plane
    normal,slip = fault_plane(t,n,p,verbose=verbose)

    # strike & dip from normal
    #
    #[strike,dip] = norm2strikedip(normal);
    #
    # make sure normals are always pointing up at hanging block
    # - this is needed to preserve the appropriate strike & dip relationship while keeping the dip in the 0-90deg range
    # - also keeps the normal vector & slip vector relationship between the fault & auxiliary planes

    # strike (in [0,360] degrees)
    v_north = normal[0]
    v_east = normal[1]
    v_up = normal[2]

    strike = atan2(-v_north,v_east)

    # for vector pointing purely in vertical direction, round-off errors might lead to wrong strike
    # (case when dip is almost zero)
    if abs(v_north) < 1.e-9 and abs(v_east) < 1.e-9:
        print("WARNING: strike is poorly constrained, fault plane normal points vertically!")
        print("         strike will be set to zero, pointing north. this will also affect the rake.\n")
        strike = 0.0

    # dip (in [0,90] degrees)
    dip = acos(v_up)

    # rake angle from slip vector and horizontal in-plane vector (in [-180,180] degree)
    v_north = slip[0]
    v_east = slip[1]
    v_up = slip[2]

    # note: acos will return values [0,pi]
    rake = acos(cos(strike) * v_north + sin(strike) * v_east)

    # return in degrees
    strike *= 180.0/pi
    dip *= 180.0/pi
    rake *= 180.0/pi

    # strike: positive clockwise from North
    # dip   : positive downward from the horizontal
    # rake  : positive counter-clockwise in the fault plane from the strike direction

    # fix sense of orientation as the above only returns positive angles
    if v_up < 0.0:
        rake = -rake

    # limits range:
    # strike in [0,360] (from North, clockwise)
    if strike < 0.0: strike += 360.0
    if strike > 360.0: strike -= 360.0

    # rather have strike = 0 than strike = 360
    if abs(strike - 360.0) < 1.e-9: strike = 0.0

    # dip in [0,90]
    if dip < 0.0: dip += 90.0
    if dip > 90.0: dip -= 90.0

    # rake in [-180,180]
    if rake < -180.0: rake += 360.0
    if rake > 180.0: rake -= 360.0

    return strike,dip,rake


#--------------------------------------------------------------------------------------------------

def convert_MT_to_SDR(mt,verbose=False):
    """
    converts moment tensor to strike-dip-rake
    """
    # decomposes the moment tensors in mt to get the maximum
    # double-couple component and then converts that to strike, dip & rake
    #
    # strike: positive clockwise from North
    # dip   : positive downward from the horizontal
    # rake  : positive counter-clockwise in the fault plane from the strike direction
    #
    # see, e.g.: https://github.com/g2e/seizmo/blob/master/cmt/mt2sdr.m

    # 1) decompose into maximum double-couple + clvd
    #    (decomposes moment tensor into isotropic + deviatoric components, double-couple and CLVD)
    mt_iso,mt_dev,mt_DC,mt_CLVD = decompose_moment_tensor(mt,verbose=verbose)

    # 2) convert to principal axes in vpa
    t,n,p = principal_axis(mt_DC,verbose=verbose)

    # 3) convert principal axes to strike/dip/rake
    strike,dip,rake = convert_principal_axis_to_SDR(t,n,p,verbose=verbose)

    return strike,dip,rake


#--------------------------------------------------------------------------------------------------


def plot_beachball(filename=None,strike=None,dip=None,rake=None,M0=None):
    """
    reads in CMT solution and plots a beachball
    """
    # obspy uses moment tensor as in
    # mt = [0.91, -0.89, -0.02, 1.78, -1.55, 0.47]

    # reads in CMT solution
    if not filename is None:
        #
        # note: obspy.read_events(filename) could also read in CMTSOLUTION format:
        #  > cat = obspy.read_events(filename,format="CMTSOLUTION")
        #  > print(cat)
        # however it assumes a correct Flinn-Engdahl region which might not hold for UTM coordinates or test CMTs
        # thus reading in manually
        moment_tensor = read_CMT(filename)

    if not strike is None:
        print("given strike/dip/rake/M0 = ",strike,dip,rake,M0)
        print("")
        # limits range:
        # strike in [0,360] (from North, clockwise)
        if strike < 0.0 or strike > 360.0:
            print("invalid strike, must be between 0 and 360 degrees")
            sys.exit(1)
        # dip in [0,90]
        if dip < 0.0 or dip > 90.0:
            print("invalid dip, must be between 0 and 90 degrees")
            sys.exit(1)
        # rake in [-180,180]
        if rake < -180.0 or rake > 180.0:
            print("invalid rake, must be between -180 and 180 degrees")
            sys.exit(1)

        # angle to rad
        strike *= pi/180.0
        dip *= pi/180.0
        rake *= pi/180.0

        # converts to moment tensor
        mt = convert_SDR_to_MT(strike,dip,rake,M0)
        moment_tensor = [mt]

    print("moment tensor MT (Mrr Mtt Mpp Mrt Mrp Mtp):")
    print(moment_tensor)
    print("")

    # moment-tensor array can hold several single CMTs.
    # we will plot for each one a beachball
    for i in range(len(moment_tensor)):
        print("moment tensor ",i+1,":")
        mt =  moment_tensor[i]

        # scalar moment M0
        M0 = get_scalar_moment(mt)

        # moment magnitude Mw
        Mw = get_moment_magnitude(mt)

        print("  magnitude of the source:")
        print("         scalar moment M0 = ",M0,"dyne-cm")
        print("      moment magnitude Mw = ",Mw)
        print("")

        # strike, dip, and rake info
        strike,dip,rake = convert_MT_to_SDR(mt,verbose=True)
        print("  strike/dip/rake = ",strike,"/",dip,"/",rake)
        print("")

        # plotting
        print("  plotting beachball")
        fig = beachball(mt, size=200, linewidth=2, facecolor='r')

        # saves figure
        if len(moment_tensor) == 1:
            outfile = "beachball" + ".png"
        else:
            outfile = "beachball_{}".format(i+1) + ".png"

        fig.savefig(outfile)
        print("  saved to: ",outfile)

    # statistics
    #for i in range(0,36+1):
    #    strike = i*10.0
    #    for j in range(0,9+1):
    #        dip = j*10.0
    #        for k in range(-18,18+1):
    #            rake = k*10.0
    #            mt = convert_SDR_to_MT(strike,dip,rake,M0)
    #            strike_out,dip_out,rake_out = convert_MT_to_SDR(mt)
    #            print("    strike/dip/rake = ",strike,"/",dip,"/",rake)
    #            print("out strike/dip/rake = ",strike_out,"/",dip_out,"/",rake_out)


    return

#--------------------------------------------------------------------------------------------------


def usage():
    print("usage: ./plot_beachball.py CMTSOLUTION or strike ")
    print("       or")
    print("       ./plot_beachball.py strike dip rake M0")
    print("   where")
    print("       CMTSOLUTION     - a source file, e.g., in DATA/CMTSOLUTION")
    print("       strike,dip,rake - source angle values (in degree, e.g. 10.0 90.0 30.0)")
    print("                         ranges: strike between [0,360], dip between [0,90] and rake between [-180,180]")
    print("       M0              - scalar moment (in dyne-cm)")

    return


if __name__ == '__main__':
    # gets argument
    if len(sys.argv) != 2 and len(sys.argv) != 5:
        usage()
        sys.exit(1)
    else:
        if len(sys.argv) == 2:
            file = sys.argv[1]
            # main routine
            plot_beachball(filename=file)
        elif len(sys.argv) == 5:
            strike = float(sys.argv[1])
            dip = float(sys.argv[2])
            rake = float(sys.argv[3])
            M0 = float(sys.argv[4])
            # main routine
            plot_beachball(strike=strike,dip=dip,rake=rake,M0=M0)



