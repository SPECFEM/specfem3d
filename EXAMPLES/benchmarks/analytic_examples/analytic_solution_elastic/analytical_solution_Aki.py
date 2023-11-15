#!/usr/bin/env python
#
#
# computes analytical solution for a double-couple source
# according to Aki & Richards, eq. 4.32 and 4.33
#
#
from __future__ import print_function

import sys

from numpy import sqrt,fabs,sin,cos,arccos,arctan2,pi,zeros
from math import erf

##############################################################################
## parameters
# source-receiver position
xs = 0.0
ys = 0.0
zs = 0.0

# receiver
xr = 6000.0
yr = 1000.0
zr = 4000.0

# (optional) file name prefix for receiver station comparison
station_prefix = "DB.Z1.FX"     # e.g., "DB.Z5.FX" -> DB.Z5.FXX.semd,DB.Z5.FXY.semd,DB.Z5.FXZ.semd
station_ending = ".semd"

# medium parameter
rho   = 2300.0
beta  = 1500.0
alpha = 2800.0

# source time function
hdur = 2.0 / 1.628    # hdur/gauss (according to SPECFEM definition)
dt   = 0.001          # time step size
M0   = 1.e16          # double-couple scalar moment (N-m)

# trace start/end time
t0 = -4.0                 # SPECFEM simulation times
t1 = (14000-1) * dt + t0  # 14000 time steps

##############################################################################


def comp_source_time_function(t,hdur):
    # quasi Heaviside, small Gaussian moment-rate tensor with hdur
    # (according to SPECFEM definition)
    val = 0.5 * (1.0 + erf(t/hdur))
    return val


def comp_source_time_function_diff(t,hdur,dt):
    # derivative of the source time function
    val_diff = (comp_source_time_function(t+dt,hdur)-comp_source_time_function(t-dt,hdur)) / (2.0*dt)
    return val_diff


def comp_source_time_function_conv(t,hdur,dt,r,alpha,beta):
    # convolution of the source time function between times r/alpha and r/beta
    nmin = int((r/alpha)/dt)
    nmax = int((r/beta)/dt)
    val_conv = 0.0
    for i in range(nmin, nmax+1):
        time_temp = i * dt
        val_conv += comp_source_time_function(t-time_temp,hdur) * time_temp * dt
    return val_conv

def analytical_solution():
    global xs,ys,zs
    global xr,yr,zr
    global station_prefix,station_ending
    global rho,alpha,beta
    global hdur,dt,M0
    global t0,t1

    # spherical coordinates
    r = sqrt((xs-xr)**2 + (ys-yr)**2 + (zs-zr)**2)
    theta = arccos((zr-zs)/r)
    phi = arctan2((yr-ys),(xr-xs))

    # user output
    print('Aki & Richards - double-couple solution')
    print('position:')
    print('  xr/yr/zr    = ',xr, yr, zr)
    print('  r/theta/phi = ',r, theta, phi)
    print('')
    print('  back-rotated x/y/z = ',r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta))
    print('')
    print('hdur:',hdur)
    print('')

    # trace start/end time
    if t1 > 0.0 and t1 > t0:
        tmin = t0
        tmax = t1
    else:
        # default: around arrivals +/- 5 * hdur
        tmin = r/alpha - 5.0 * hdur
        tmax = r/beta + 5.0 * hdur
    ntmin = int(tmin/dt)
    ntmax = int(tmax/dt)

    print('trace length:')
    print('  tmin / ntmin :',tmin, ntmin)
    print('  tmax / ntmax :',tmax, ntmax)
    print('')

    #compute factors in the AKI & RICHARDS solution
    cn  = 1.0/(4.0 * pi * rho * r**4) * M0
    cip = 1.0/(4.0 * pi * rho * alpha**2 * r**2) * M0
    cis = 1.0/(4.0 * pi * rho * beta**2 * r**2) * M0
    cfp = 1.0/(4.0 * pi * rho * alpha**3 * r) * M0
    cfs = 1.0/(4.0 * pi * rho * beta**3 * r) * M0

    # Aki & Richards, eqs. (4.33)
    # components index: 0 == r, 1 == theta, 2 == phi
    an = zeros(3)
    an[0] = 9.0 * sin(2.0*theta) * cos(phi)
    an[1] = -6.0 * cos(2.0*theta) * cos(phi)
    an[2] = 6.0 * cos(theta) * sin(phi)

    aip = zeros(3)
    aip[0] = 4.0*sin(2.0*theta)*cos(phi)
    aip[1] = -2.0*cos(2.0*theta)*cos(phi)
    aip[2] = 2.0*cos(theta)*sin(phi)

    ais = zeros(3)
    ais[0] = -3.0*sin(2.0*theta)*cos(phi)
    ais[1] = 3.0*cos(2.0*theta)*cos(phi)
    ais[2] = -3.0*cos(theta)*sin(phi)

    afp = zeros(3)
    afp[0] = sin(2.0*theta)*cos(phi)
    afp[1] = 0.0
    afp[2] = 0.0

    afs = zeros(3)
    afs[0] = 0.0
    afs[1] = cos(2.0*theta)*cos(phi)
    afs[2] = -cos(theta)*sin(phi)

    usph = zeros(3)
    ucar = zeros(3)

    tmp_nf = zeros(3)
    tmp_if = zeros(3)
    tmp_ff = zeros(3)

    # displacement
    if len(station_prefix) > 0:
        name1 = station_prefix + "X" + station_ending
        name2 = station_prefix + "Y" + station_ending
        name3 = station_prefix + "Z" + station_ending
    else:
        # default naming
        name1 = "Ux.dat"
        name2 = "Uy.dat"
        name3 = "Uz.dat"

    print('trace names:')
    print('  ',name1,name2,name3)
    print('')

    # file output
    f1 = open(name1,'w')
    f2 = open(name2,'w')
    f3 = open(name3,'w')

    # format: #time #u_cartesian #u_spherical
    info = "# Aki and Richards - double-couple solution eq.(4.32) and (4.33)\n"
    info += "#\n"
    info += "# homogeneous elastic medium: rho/vp/vs   = {} / {} / {}\n".format(rho,alpha,beta)
    info += "# receiver station          : x/y/z       = {} / {} / {}\n".format(xr,yr,zr)
    info += "#                           : r/theta/phi = {} / {} / {}\n".format(r,theta,phi)
    info += "#\n"

    comp1 = "# solution component        : cartesian X / spherical R\n"
    comp2 = "# solution component        : cartesian Y / spherical THETA\n"
    comp3 = "# solution component        : cartesian Z / spherical PHI\n"

    format = "#\n"
    format += "# format:\n"
    format += "#time \t#u_cartesian \t#u_spherical\n"

    f1.write(info + comp1 + format)
    f2.write(info + comp2 + format)
    f3.write(info + comp3 + format)

    # source time functions
    stf1 = open('stf_step_time_source','w')
    stf2 = open('stf_diff_step_time_source','w')
    stf3 = open('stf_conv_step_time_source','w')

    # format: #time #stf
    stf1.write('#time #stf\n')
    stf2.write('#time #stf_diff\n')
    stf3.write('#time #stf_conv\n')

    print('')
    print('writing wavefield solution...')

    for it in  range(ntmin, ntmax+1):
        time = it * dt

        # source time functions
        stf_p = comp_source_time_function(time-r/alpha,hdur)
        stf_s = comp_source_time_function(time-r/beta,hdur)

        stf_p_diff = comp_source_time_function_diff(time-r/alpha,hdur,dt)
        stf_s_diff = comp_source_time_function_diff(time-r/beta,hdur,dt)

        stf_conv = comp_source_time_function_conv(time,hdur,dt,r,alpha,beta)

        stf1.write("{} \t{} \t{}\n".format(time, stf_p, stf_s))
        stf2.write("{} \t{} \t{}\n".format(time, stf_p_diff, stf_s_diff))
        stf3.write("{} \t{}\n".format(time, stf_conv))

        # wavefield terms
        # components (r,theta,phi)
        for i in range(3):
            # near-field
            tmp_nf[i] = cn  * an[i]  * stf_conv
            # intermediate-field
            tmp_if[i] =  cip * aip[i] * stf_p + cis * ais[i] * stf_s
            # far-field
            tmp_ff[i] = cfp * afp[i] * stf_p_diff + cfs * afs[i] * stf_s_diff
            # spherical solution
            usph[i] = tmp_nf[i] + tmp_if[i] + tmp_ff[i]

        ## displacement vector in spherical coordinates d = v_r * unit_r + v_t * unit_theta + v_p * unit_phi
        #v_r = usph(1)
        #v_t = usph(2)
        #v_p = usph(3)
        ## rotation to cartesian system
        ## see: https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates
        ## displacement vector in cartesian coordinates d = v_x * unit_x + v_y * unit_y + v_z * unit_z
        #v_x = v_r * sin(theta)*cos(phi) + v_t * cos(theta)*cos(phi) + v_p * (-sin(phi))
        #v_y = v_r * sin(theta)*sin(phi) + v_t * cos(theta)*sin(phi) + v_p * cos(phi)
        #v_z = v_r * cos(theta)          + v_t * (-sin(theta))      !+ v_p * 0.0

        ## coordinate system convention:
        # SPECFEM:
        # the x axis points East
        # the y axis points North
        # the z axis points up
        #
        # Aki & Richards:
        # the x axis points North
        # the y axis points East
        # the z axis points down
        #
        # results for:
        # unit_r     = (sin(theta) cos(phi), sin(theta) sin(phi), cos(theta)  )
        # unit_theta = (cos(theta) cos(phi), cos(theta) sin(phi), - sin(theta))
        # unit_phi   = (-sin(phi)          , cos(phi)           , 0           )
        #
        # slip in (x1,x2) plane along x1
        #
        ## Aki & Richards convention:
        #v_x = v_r * sin(theta)*cos(phi) + v_t * cos(theta)*cos(phi) + v_p * (-sin(phi))
        #v_y = v_r * sin(theta)*sin(phi) + v_t * cos(theta)*sin(phi) + v_p * cos(phi)
        #v_z = v_r * cos(theta)          + v_t * (-sin(theta))      !+ v_p * 0.0
        #
        ## conversion to SPECFEM?
        #
        # STATIONS position:
        # #name   #network #y      #x      #(ignored) # z(Par_file: USE_SOURCES_RECEIVERS_Z = .true.)
        # X1      DB       7000.0  8000.0  0.0        4000.0
        #
        # CMTSOLUTION:
        # r -> z, theta -> -y, phi -> x
        #
        #  Mxx = Mpp
        #  Myy = Mtt
        #  Mzz = Mrr
        #  Myz = -Mrt
        #  Mxz = Mrp
        #  Mxy = -Mtp
        #
        # the setup in Aki & Richards corresponds to component Mrp being non-zero
        # CMTSOLUTIONs use dyne-cm, where as here M0 is in N-m -> conversion factor M_(dyne-cm) = 10**7 M_(N-m)
        #
        ## convert vector field in r/theta/phi to cartesian coordinate system x/y/z
        vec_r = usph[0]
        vec_t = usph[1]
        vec_p = usph[2]

        # u(r,theta,phi) -> u(x,y,z)
        ucar[0] = vec_r * sin(theta)*cos(phi) + vec_t * cos(theta)*cos(phi) - vec_p * sin(phi)
        ucar[1] = vec_r * sin(theta)*sin(phi) + vec_t * cos(theta)*sin(phi) + vec_p * cos(phi)
        ucar[2] = vec_r * cos(theta)          - vec_t * sin(theta)

        # format: #time #u_cartesian #u_spherical
        f1.write("{} \t{} \t{} \t{} {} {}\n".format(time,ucar[0],usph[0],tmp_nf[0],tmp_if[0],tmp_ff[0]))   # #t #x #r
        f2.write("{} \t{} \t{} \t{} {} {}\n".format(time,ucar[1],usph[1],tmp_nf[1],tmp_if[1],tmp_ff[1]))   # #t #y #theta
        f3.write("{} \t{} \t{} \t{} {} {}\n".format(time,ucar[2],usph[2],tmp_nf[2],tmp_if[2],tmp_ff[2]))   # #t #z #phi

    f1.close()
    f2.close()
    f3.close()

    stf1.close()
    stf2.close()
    stf3.close()

    print('')
    print('written to: ')
    print('  ',name1)
    print('  ',name2)
    print('  ',name3)
    print('')
    print('done')


def usage():
    print('usage: ./analytical_solution.py xr yr zr station_prefix')
    print('  with')
    print('     xr,yr,zr       - receiver position in m (e.g., 8000.0 1000.0 4000.0)')
    print('     station_prefix - station name prefix (e.g., DB.Z1.FX to compare with SPECFEM)')
    print('')


if __name__ == '__main__':
    # gets (optional) arguments
    if len(sys.argv) > 1:
        if len(sys.argv) != 5:
            usage()
            sys.exit(1)

        # receiver position
        xr = float(sys.argv[1])
        yr = float(sys.argv[2])
        zr = float(sys.argv[3])

        # station prefix
        station_prefix = sys.argv[4]

    # solution
    analytical_solution()
