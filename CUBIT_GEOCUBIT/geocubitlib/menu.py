#############################################################################
# menu.py
# this file is part of GEOCUBIT                                             #
#                                                                           #
# Created by Emanuele Casarotti                                             #
# Copyright (c) 2008 Istituto Nazionale di Geofisica e Vulcanologia         #
#                                                                           #
#############################################################################
#                                                                           #
# This program is free software; you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# This program is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License along   #
# with this program; if not, write to the Free Software Foundation, Inc., #
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.               #
#                                                                           #
#############################################################################
from __future__ import print_function

import getopt
import sys

from utilities import get_cubit_version


def usage():
    txt = """
GEOCUBIT HELP...

1) UTILITIES

 check the configuration of the libraries and dependencies:
 GEOCUBIT.py --chklib

 check the parameter file:
 GEOCUBIT.py --chkcfg --cfg = [filename]

2) CREATE GEOMETRY

 create a surface for regular ascii grid or ascii lines defining a skin
 GEOCUBIT.py --surface = [surface file] (--regulargrid = [...] --skin = [...])

 create a plane surface
 GEOCUBIT.py --plane --x1 = [x, y, z]
                     --x2 = [x, y, z]
                     --x3 = [x, y, z]
                     --x4 = [x, y, z]
                     --unit = [utm/geo]

 create acis surfaces using a parameter file:
 GEOCUBIT.py --build_surface --cfg = [filename]

 create cubit volumes using a parameter file:
 GEOCUBIT.py --build_volume --cfg = [filename] (--id_proc = [num_processor])


3) MESHING

meshing a volumes
 GEOCUBIT.py --mesh --cfg = [filename] (--id_proc = [num_processor])
 - note: without the --build_volume flag the script recall
 an old 'geometry_vol_[id_proc].cub' file

 build a volume and mesh it....
 GEOCUBIT.py --build_volume --mesh --cfg = [filename]
             (--id_proc = [num_processor, default = 0])

4) FINALIZING AND EXPORTING

 collect some cubit files and merge in a single free mesh cubitfile
 GEOCUBIT.py --collect --merge --meshfiles = [files] --cpux = N --cpuy = N
             (--rangecpux = [cpuxmin, cpuxmax],
             --rangecpuy = [cpuymin, cpuymax]
             --output = [YourMeshName]
             --outdir = [YourDir]
             --step_tolerance=1
             --save_cubfile
             --starting_tolerance = 100)


 export a single cubit mesh file (with defined blocks) in a SPECFEM3D mesh
 GEOCUBIT.py --export2SPECFEM3D --meshfiles = [filename]
             (--listblock = block1, block2, .., blockN
              --listflag = [specfem flag, i.e. --listflag = 1, 2, 3, -1])
              --SEMoutput = [YourOutputDir])

"""
    print(txt)


# print('reading options....')
try:
    if hasattr(sys, 'argv'):
        opts, args = getopt.getopt(sys.argv[1:], "sjmohbp1",
                                   ["starting_tolerance=", "save_cubfile",
                                    "step_tolerance=", "hex27", "cpml_size=",
                                    "top_absorbing", "cpml", "decimate",
                                    "addsea", "SEMoutput=", "qlog", "mfast",
                                    "curverefining=", "output=", "rangecpux=",
                                    "rangecpuy=", "equivalence", "listflag=",
                                    "listblock=", "cpux=", "cpuy=",
                                    "exofiles=", "partitioner", "plane",
                                    "x1 = ", "x2 = ", "x3 = ", "x4 = ",
                                    "unit=", "chkcfg", "mat=",
                                    "merge_tolerance=", "export2SPECFEM3D",
                                    "mesh", "chklib", "cfg=", "job = ",
                                    "basin", "help", "id_proc=", "surface=",
                                    "script", "jou", "strat", "MPI",
                                    "regulargrid=", 'skin=', "build_surface",
                                    "build_volume", "merge1", "merge2",
                                    "merge", "collect", "meshfiles="])

except Exception as e:
    print(e)
    sys.exit()

output = 'totalmesh_merged'
SPECFEM3D_output_dir = '.'
verbose = False
surface = False
script = False
mesh = False
basin = False
configuration = None
id_proc = 0
single = False
nomesh = False
jobid = 0
build_surface = False
build_volume = False
meshing = False
ckbound_method1 = False
ckbound_method2 = False
collect = False
export2SPECFEM3D = False
merge_tolerance = 0
starting_tolerance = 100
material_file = 'material_archive.dat'
material_assignement = []
chkcfg = False
create_plane = False
create_partitioner = False
cubfiles = None
exofiles = None
listflag = None
listblock = None
cpux = 1
cpuy = 1
cpuxmin = 0
cpuymin = 0
cpuxmax = None
cpuymax = None
curverefining = False
add_sea = False
decimate = False

cpml = False
cpml_size = False
top_absorbing = False

qlog = False
hex27 = False

check_merging = False
save_cubfile = False

cubit_version = get_cubit_version()


if opts:
    for o, value in opts:
        if o in ('--starting_tolerance'):
            starting_tolerance = float(value)
        # print(o, value)
        if o in ('--save_cubfile'):
            save_cubfile = True
        if o in ('--step_tolerance'):
            step_tolerance = float(value)
        if o in ('--hex27'):
            hex27 = True

        if o in ('--cpml'):
            cpml = True
            if '--cpml_size' in o:
                cpml = True
                cpml_size = float(value)
        if o in ('--top_absorbing'):
            cpml = True
            top_absorbing = True
        if '--cpml_size' in o:
            cpml = True
            cpml_size = float(value)
        if o in ('--decimate'):
            decimate = True
        if o in ('--partitioner'):
            create_partitioner = True
        if o == ('--surface'):
            surface = True
            surface_name = value
        if o == ('--build_surface'):
            build_surface = True
        if o == ('--build_volume'):
            build_volume = True
        if surface and o == ('--regular_grid'):
            surface_type = 'regular_grid'
            tmp = value.split('/')
            if len(tmp) == 4:
                num_x, num_y, unit, delimiter = value.split('/')
            elif len(tmp) == 3:
                num_x, num_y, unit = value.split('/')
            num_x = int(num_x)
            num_y = int(num_y)
            delimiter = ' '
        if surface and o == ('--skin'):
            surface_type = 'skin'
            tmp = value.split('/')
            if len(tmp) == 4:
                directionx, directiony, unit, delimiter = value.split('/')
            elif len(tmp) == 3:
                directionx, directiony, unit = value.split('/')
            directiony = int(directiony)
            directionx = int(directionx)
        if o in ('--plane'):
            create_plane = True
            cfg_name = False
        if o in ('--x1'):
            x1 = value
        if o in ('--x2'):
            x2 = value
        if o in ('--x3'):
            x3 = value
        if o in ('--x4'):
            x4 = value
        if o in ('--unit'):
            unit = value
        #
        if o in ('--build_volume'):
            build_volume = True
        if o in ("-h", "--help"):
            usage()
            sys.exit(2)
        # check the configuration
        if o in ("--chklib"):
            import start as start
            mpiflag, iproc, numproc, mpi = start.start_mpi()
            if mpiflag:
                print('--------, MPI ON, parallel mesher ready')
            else:
                print('--------, MPI OFF, serial mesher ready')
            numpy = start.start_numpy()
            print('--------, Numpy ON')
            cubit = start.start_cubit()
            print('--------, CUBIT ON')
            sys.exit()
        if o in ("--cfg"):
            cfg_name = value
            try:
                if open(cfg_name):
                    pass
            except IOError as e:
                print('error opening ', cfg_name)
                print(e)
                import sys
                sys.exit()
        if o == ('--surface'):
            surface = True
            surface_name = value
        if o == ("--mesh"):
            meshing = True
        if o in ("-1"):
            single = True
        if o in ("--collect"):
            collect = True
        if o in ("--merge2") and o != '--merge':
            ckbound_method2 = True
        if o in ("--equivalence", "--merge1", "--merge"):
            ckbound_method1 = True
        if o in ("--mfast"):
            ckbound_method1 = True
            ckbound_method2 = True
        if o == ("--meshfiles"):
            cubfiles = value
            import glob
            nf = glob.glob(cubfiles)
            if len(nf) > 0:
                print('cubfiles ', nf)
            else:
                print('files not found: ', cubfiles)
                import sys
                sys.exit()
        if o in ("--exofiles"):
            exofiles = value
        if o in ("--export2SPECFEM3D"):
            export2SPECFEM3D = True
        if o in ("--merge_tolerance"):
            if o != '--merge' and o != '--merge2' and o != '--merge1':
                merge_tolerance = map(float, value.split(','))
        if o in ("--mat"):
            material_assignement.append(
                [value.split(',')[0], value.split(',')[1]])
        if o in ("--listblock"):
            listblock = map(int, value.split(','))
        if o in ("--listflag"):
            listflag = map(int, value.split(','))
        if o in ("--chkcfg"):
            chkcfg = True
        if o in ('--id_proc'):
            id_proc = int(value)
        if o in ('--rangecpux'):
            cpuxmin = int(value.split(',')[0])
            if cubit_version >= 14.0:
                cpuxmax = int(value.split(',')[1])
            else:
                cpuxmax = int(value.split(',')[1]) + 1
        if o in ('--rangecpuy'):
            if cubit_version >= 14.0:
                cpuymax = int(value.split(',')[1])
            else:
                cpuymax = int(value.split(',')[1]) + 1
            cpuymin = int(value.split(',')[0])
        if o in ('--cpux'):
            cpux = int(value)
        if o in ('--cpuy'):
            cpuy = int(value)
        if o in ("--qlog"):
            qlog = True
        if o in ("--output", "-o"):
            output = value
        if o in ("--curverefining"):
            curverefining = value.split(',')
        if o in ("--SEMoutput"):
            SPECFEM3D_output_dir = value
            import os
            try:
                os.makedirs(SPECFEM3D_output_dir)
            except OSError:
                pass
        if o in ("--addsea"):
            add_sea = True
    print(cpuxmax, cpuymax)
    if cpuymax:
        pass
    elif cpuy > 1:
        if cubit_version >= 14.0:
            cpuymax = cpuy - 1
        else:
            cpuymax = cpuy
    else:
        cpuymax = 1
    if cpuxmax:
        pass
    elif cpux > 1:
        if cubit_version >= 14.0:
            cpuxmax = cpux - 1
        else:
            cpuxmax = cpux
    else:
        cpuxmax = 1
    print(cpuxmax, cpuymax)

    if cpml:
        if not cpml_size:
            print('specify the size of the cpml boundaries')
            import sys
            sys.exit()
        elif cpml_size <= 0:
            print('cpml size negative, please check the parameters')
            import sys
            sys.exit()

    if chkcfg is True:
        import start as start
        cfg = start.start_cfg()
        d = cfg.__dict__
        ks = d.keys()
        ks.sort()
        for k in ks:
            if '__' not in k and '<' not in str(d[k]) and d[k] is not None:
                txt = str(k) + ' -----> ' + \
                    str(d[k])
                txt = txt.replace("'", "").replace('"', '')
                print(txt)
    else:
        try:
            import start as start
            cfg = start.start_cfg()
            f = open('cfg.log', 'w')
            print>>f, 'CFG FILE: ', cfg_name
            d = cfg.__dict__
            ks = d.keys()
            ks.sort()
            for k in ks:
                if '__' not in k and '<' not in str(d[k]) and d[k] is not None:
                    txt = str(k) + ' -----> ' + str(d[k])
                    txt = txt.replace("'", "").replace('"', '')
                    print>>f, txt
            f.close()
        except:
            pass
elif opts == []:
    print(__name__)
    usage()
    import sys
    sys.exit()
