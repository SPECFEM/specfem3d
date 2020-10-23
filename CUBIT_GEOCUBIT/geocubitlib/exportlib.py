#!/usr/bin/env python
#############################################################################
# exportlib.py
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
# This program is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License along   #
# with this program; if not, write to the Free Software Foundation, Inc.,   #
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.               #
#                                                                           #
#############################################################################
from __future__ import print_function

try:
    import start as start
    cubit = start.start_cubit()
except:
    try:
        import cubit
    except:
        print("error importing cubit, check if cubit is installed")
        pass

import glob
from utilities import get_cubit_version
print('version 2.2')


class VersionException(Exception):
    pass


class MergingError(Exception):
    pass


def invert_dict(d):
    inv = {}
    for k, v in d.iteritems():
        keys = inv.setdefault(v, [])
        keys.append(k)
    return inv


def importing_cubfiles(cubfiles):
    import re
    rule_st = re.compile("(.+)_[0-9]+\.")
    rule_ex = re.compile(".+_[0-9]+\.(.+)")
    rule_int = re.compile(".+_([0-9]+)\.")
    filenames = glob.glob(cubfiles)
    try:
        st = rule_st.findall(filenames[0])[0]
        ex = rule_ex.findall(filenames[0])[0]
        listflag = True
    except:
        ex = ''
        listflag = False
    if ex == 'cub':
        cubflag = True
    else:
        cubflag = False
    list_int = []
    fs = []
    try:
        for f in filenames:
            i = int(rule_int.findall(f)[0])
            list_int.append(i)
        list_int.sort()
        for i, ind in enumerate(list_int):
            f = st + '_' + str(ind) + '.' + ex
            fs.append(f)
    except:
        pass
    if listflag:
        filenames = fs
    else:
        pass
    return len(filenames), list_int, filenames, cubflag


def collect_new(cpuxmin=0, cpuxmax=1, cpuymin=0, cpuymax=1, cpux=1, cpuy=1,
                cubfiles=False, ckbound_method1=False, ckbound_method2=False,
                merge_tolerance=None, curverefining=False,
                outfilename='totalmesh_merged', qlog=False,
                export2SPECFEM3D=False, listblock=None, listflag=None,
                outdir='.', add_sea=False, decimate=False, cpml=False,
                cpml_size=False, top_absorbing=False, hex27=False,
                save_cubfile=False, check_merging=False,
                starting_tolerance=100.):
    #
    cubit.cmd('set info off')
    cubit.cmd('set echo off')
    cubit.cmd('set journal off')
    # cubit.cmd('set error off')
    version_cubit = get_cubit_version()

    print('##########################################################')
    print('#')
    print('# collecting mesh')
    print('#')
    print('##########################################################')

    if decimate:
        if version_cubit >= 14.0:
            raise VersionException('check cubit version, decimate capability \
                                   has been tested only with cubit <= 12.2')

    if version_cubit <= 12.2:
        collecting_merging(cpuxmin, cpuxmax, cpuymin, cpuymax, cpux, cpuy,
                           cubfiles=cubfiles, ckbound_method1=ckbound_method1,
                           ckbound_method2=ckbound_method2,
                           merge_tolerance=merge_tolerance, decimate=decimate)
    elif version_cubit >= 14:
        collecting_merging_new(cpuxmin, cpuxmax + 1, cpuymin, cpuymax + 1,
                               cpux, cpuy, cubfiles=cubfiles,
                               check_merging=check_merging,
                               starting_tolerance=starting_tolerance)
    else:
        raise VersionException('check cubit version, parallel capability \
                               of geocubit is working with \
                               cubit/trelis 14 or later (or cubit 12.2)')

    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    cubit.cmd('set journal on')
    # cubit.cmd('set error on')

    ##
    ## refinements
    ##
    if curverefining:
        if version_cubit <= 12.2:
            block = 1001  # topography
            refine_closecurve(block, curverefining, acis=True)
        else:
            raise VersionException(
                'check cubit version, refine curve capability has been tested \
                only with cubit <= 12.2')

    ##
    ## sea
    ##
    if add_sea:
        if version_cubit <= 12.2:
            block = 1001
            add_sea_layer(block=block)
        else:
            raise VersionException(
                'check cubit version, sea capability has been tested \
                only with cubit <= 12.2')

    # output directory
    outdir2 = '/'.join(x for x in outfilename.split('/')[:-1])
    if outdir2 == '':
        outdir2 = outdir + '/'
    else:
        outdir2 = outdir + '/' + outdir2 + '/'

    import os
    try:
        os.makedirs(outdir2)
    except OSError:
        pass

    # mesh file output
    cubit.cmd('compress all')
    command = "export mesh '" + outdir2 + outfilename + \
        ".e' block all overwrite xml '" + outdir2 + outfilename + ".xml'"
    cubit.cmd(command)

    # outputs block infos
    f = open(outdir2 + 'blocks.dat', 'w')
    blocks = cubit.get_block_id_list()
    #
    for block in blocks:
        name = cubit.get_exodus_entity_name('block', block)
        element_count = cubit.get_exodus_element_count(block, "block")
        nattrib = cubit.get_block_attribute_count(block)
        attr = [cubit.get_block_attribute_value(
            block, x) for x in range(0, nattrib)]
        ty = cubit.get_block_element_type(block)
        f.write(str(block) + '; ' + name + ' ; nattr ' + str(nattrib) + ' ; ' +
                ' '.join(str(x) for x in attr) + ' ; ' + ty + ' ' +
                str(element_count) + '\n')
    f.close()

    #
    #
    cubit.cmd('set info echo journ off')
    cmd = 'del group all'
    cubit.silent_cmd(cmd)
    cubit.cmd('set info echo journ on')
    #
    #
    print('end meshing')
    #
    #
    if qlog:
        print('\n\nQUALITY CHECK.... ***************\n\n')
        import quality_log
        tq = open(outdir2 + outfilename + '.quality', 'w')
        max_skewness, min_length = quality_log.quality_log(tq)
    #
    #
    #
    if export2SPECFEM3D:
        e2SEM(files=False, listblock=listblock,
              listflag=listflag, outdir=outdir,
              cpml=cpml, cpml_size=cpml_size,
              top_absorbing=top_absorbing, hex27=hex27)

    if save_cubfile:
        vol_blocks = [x for x in blocks if x <= 1000]
        cubit.cmd("create mesh geometry block " +
                  ' '.join(str(x) for x in vol_blocks) +
                  " feature_angle 135.0")
        command = "save as '" + outdir2 + outfilename + ".cub' overwrite"
        print(command)
        cubit.cmd(command)


def e2SEM(files=False, listblock=None, listflag=None, outdir='.',
          cpml=False, cpml_size=False, top_absorbing=False, hex27=False):
    import glob
    if files:
        filenames = glob.glob(files)
        for f in filenames:
            print(f)
            extension = f.split('.')[-1]
            if extension == 'cub':
                cubit.cmd('open "' + f + '"')
            elif extension == 'e':
                cubit.cmd('import mesh "' + f + '" no_geom')
            else:
                print(extension)
    if listblock and listflag:
        pass
    else:
        listblock = []
        listflag = []
        block_list = list(cubit.get_block_id_list())
        for block in block_list:
            ty = cubit.get_block_element_type(block)
            if 'HEX' in ty:
                listblock.append(block)
                # listflag.append(block)
        listflag = range(1, len(block_list) + 1)
    #
    for ib, iflag in zip(listblock, listflag):
        cubit.cmd("block " + str(ib) + " attribute count 1")
        cubit.cmd("block " + str(ib) + " attribute index 1 " + str(iflag))
    #
    import cubit2specfem3d
    cubit2specfem3d.export2SPECFEM3D(outdir, cpml=cpml, cpml_size=cpml_size,
                                     top_absorbing=top_absorbing, hex27=hex27)


def collecting_block(store_group_name, iproc=0, xmin=[0], xmax=[0], ymin=[0],
                     ymax=[0], index_block=0):
    block_list = list(cubit.get_block_id_list())
    block_list.sort()
    block_hex = [x for x in block_list if x <= 1000]
    block_side = [x for x in block_list if x > 1000]
    #print('collecting block list = ',block_list)
    #print('  hex list  = ',block_hex)
    #print('  side list = ',block_side)
    print('collecting block: iproc = ',iproc,xmin,xmax,ymin,ymax)
    print('  block_hex = ',block_hex,'block_side = ',block_side)
    ## hexahedrals
    for ib, block in enumerate(block_hex):
        #print('  hex ib = ',ib,'block = ',block)
        if index_block == 0:
            cubit.cmd("group 'vol" + str(block) +
                      "' add Hex in block " + str(block))
            store_group_name.append('vol' + str(block))
            cubit.cmd("del block " + str(block))
        else:
            cubit.cmd("group '" + store_group_name[ib] +
                      "' add Hex in block " + str(block))
            cubit.cmd("del block " + str(block))
    ## faces
    for ib, side in enumerate(block_side):
        #print('  faces ib = ',ib,'side = ',side)
        if side == 1004:
            if iproc in ymin:
                cubit.cmd("group 'ymin' add face in block " + str(side))
            else:
                cubit.cmd("group 'lateral' add face in block " + str(side))
        elif side == 1003:
            if iproc in xmin:
                cubit.cmd("group 'xmin' add face in block " + str(side))
            else:
                cubit.cmd("group 'lateral' add face in block " + str(side))
        elif side == 1006:
            if iproc in ymax:
                cubit.cmd("group 'ymax' add face in block " + str(side))
            else:
                cubit.cmd("group 'lateral' add face in block " + str(side))
        elif side == 1005:
            if iproc in xmax:
                cubit.cmd("group 'xmax' add face in block " + str(side))
            else:
                cubit.cmd("group 'lateral' add face in block " + str(side))
        elif side == 1001:
            cubit.cmd("group 'topo' add face in block " + str(side))
        elif side == 1002:
            cubit.cmd("group 'bot' add face in block " + str(side))
        cubit.cmd("del block " + str(side))

    ## check if faces in group lateral
    ilateral = cubit.get_id_from_name('lateral')
    lateral_nodes = cubit.get_group_nodes(ilateral)
    print('  lateral nodes: ', len(lateral_nodes))

    return store_group_name


def check_lateral_nodes(name_group='lateral'):
    cubit.cmd("group 'lateral_nodes' add Node in face in group " + name_group)
    ilateral_nodes = cubit.get_id_from_name('lateral_nodes')
    lateral_nodes = cubit.get_group_nodes(ilateral_nodes)
    cubit.cmd('del group ' + str(ilateral_nodes))
    print(name_group, ' nodes ', len(lateral_nodes))
    return lateral_nodes


def prepare_equivalence_new(name_group='lateral'):
    print('equivalence group ',name_group)
    length = {}
    cmd = "group 'tmpn' add edge in face in group " + name_group
    cubit.cmd(cmd)
    ge = cubit.get_id_from_name("tmpn")
    e1 = cubit.get_group_edges(ge)
    lengthmin = 1e9
    for e in e1:
        lengthmin = min(lengthmin, cubit.get_mesh_edge_length(e))
        length[e] = lengthmin * 0.5
    cubit.cmd('delete group ' + str(ge))

    #print('  equivalence edge lengths ',length)
    if len(length) > 0:
        minvalue = min(length.values())
        maxvalue = max(length.values())
    else:
        minvalue = 100.
        maxvalue = 100.
    #
    print('  min length: ', minvalue, 'max length: ', maxvalue)
    if minvalue != 0:
        nbin = int((maxvalue / minvalue)) + 1
        factor = (maxvalue - minvalue) / nbin
    else:
        nbin = 0
        factor = 0.0
    dic_new = {}
    for k in length.keys():
        if factor != 0.0:
            dic_new[k] = int((length[k] - minvalue) / factor)
        else:
            dic_new[k] = 0.0
    inv_length = invert_dict(dic_new)

    ks = inv_length.keys()
    ks.sort()
    for k in range(0, len(inv_length.keys()) - 1):
        inv_length[ks[k]] = inv_length[ks[k]] + inv_length[ks[k + 1]]

    print('  edge lengths ',inv_length.keys(), factor, minvalue)

    return factor, minvalue, maxvalue, inv_length


def merging_node_new(tol, clean=True, graphic_debug=False):
    empty = False
    print('tolerance ', tol)
    cubit.cmd("topology check coincident node node in \
              group coincident_lateral_nodes tolerance " + str(tol) + " highlight brief result \
              group 'merging_lateral_nodes'")
    group_exist = cubit.get_id_from_name("merging_lateral_nodes")
    if not group_exist:
        print('no nodes in this tolerance range')
    else:
        merging_nodes = cubit.get_group_nodes(group_exist)
        if graphic_debug:
            cubit.cmd('draw group lateral')
            cubit.cmd('high group merging_lateral_nodes')
        print('merging ', len(merging_nodes), ' nodes.....')
        cubit.cmd("equivalence node in merging_lateral_nodes \
                  tolerance " + str(tol * 2))
        if clean:
            cubit.cmd("group coincident_lateral_nodes \
                      remove node in group merging_lateral_nodes")
            cubit.cmd("delete Group merging_lateral_nodes")
        ic_nodes = cubit.get_id_from_name('coincident_lateral_nodes')
        c_nodes = cubit.get_group_nodes(ic_nodes)
        print(len(c_nodes))
        if len(c_nodes) == 0:
            empty = True
        if graphic_debug:
            cubit.cmd('draw group lateral')
            cubit.cmd('high group coincident_lateral_nodes')
            cubit.cmd('quality hex all jacobian \
                      global high 0 draw mesh draw add')
        return empty


def graphic_merging(tol, step_tol=None, maxtol=None):
    """
    routine for merging chunks in cubit/trelis GUI
    tol :: tolerance starting value
    step_tol :: the value that iteratively increases tol
    maxtol :: max value of tolerance
    """
    if not step_tol:
        step_tol = tol / 10.
    if not maxtol:
        maxtol = tol * 100

    cubit.cmd('group \'coincident_lateral_nodes\' add \
              Node in face in group lateral')
    isempty = False
    while isempty:
        isempty = merging_node_new(tol, clean=True, graphic_debug=True)
        tol = tol + step_tol
        if tol > maxtol:
            print('tolerance greater than the max length of the edges, \
                  please check the mesh')


def collecting_merging_new(cpuxmin=0, cpuxmax=0, cpuymin=0, cpuymax=0, cpux=1,
                           cpuy=1, cubfiles=False, check_merging=False,
                           starting_tolerance=None, step_tolerance=None):
    # import glob
    # import re
    #
    ##
    try:
        from boundary_definition import check_bc, map_boundary
    except:
        pass
    #
    number_of_chunks = cpux * cpuy
    print('number of chunks: ', number_of_chunks)

    xmin, xmax, ymin, ymax, listfull = map_boundary(cpuxmin, cpuxmax, cpuymin, cpuymax, cpux, cpuy)
    print('list of processes with boundary:')
    print('  xmin: ', xmin)
    print('  xmax: ', xmax)
    print('  ymin: ', ymin)
    print('  ymax: ', ymax)
    print('  full list: ', listfull)
    if 1 < number_of_chunks < max(listfull):
        raise MergingError('error mapping the chunks')
    #
    if cubfiles:
        nf, listiproc, filenames, cubflag = importing_cubfiles(cubfiles)
        print(nf, listiproc, filenames, cubflag, listfull)
    else:
        nf = 0
        filenames = []
        iproc = 0
    print('cubfiles       : ',cubfiles)
    print('number of files: ',nf)
    #

    index_block = -1
    store_group_name = []
    side_name = ['topo', 'xmin', 'ymin',
                 'xmax', 'ymax', 'bot']
    side_val = ['1001', '1003', '1004',
                '1005', '1006', '1002']
    side_block_name = ['face_topo', 'face_abs_xmin', 'face_abs_ymin',
                       'face_abs_xmax', 'face_abs_ymax', 'face_abs_bottom']
    cubit.cmd('set duplicate block elements on')

    if nf > 0:
        for iproc, filename in zip(listiproc, filenames):
            print(iproc, filename, iproc in listfull)
            try:
                if iproc in listfull:
                    print(filename)
                    index_block = index_block + 1
                    if cubflag:
                        cubit.cmd('import cubit "' + filename + '"')
                    else:
                        cubit.cmd('import mesh "' + filename + '" block all  no_geom')
            except:
                cubit.cmd('import mesh "' + filename + '" block all  no_geom')
            # print(iproc,xmin,xmax,ymin,ymax,iproc in xmin,iproc in xmax,iproc in ymin,iproc in ymax)
            print("iproc ",iproc," mesh imported")

            store_tmp = collecting_block(store_group_name, iproc, xmin, xmax, ymin, ymax, index_block)
            print('store tmp ',store_tmp)
            if len(store_tmp) != 0:
                store_group_name = store_tmp
            # lateral_nodes = check_lateral_nodes()
            print('store group name ',store_group_name)

        cubit.cmd('save as "tmp_nomerging.cub" overwrite ')

    else:
        if number_of_chunks == 1:
            from geocubitlib import boundary_definition
            boundary_definition.define_bc()
        else:
            check_bc(iproc, xmin, xmax, ymin, ymax, cpux, cpuy,
                     cpuxmin, cpuxmax + 1, cpuymin, cpuymax + 1)
        cubit.cmd('disassociate mesh from volume all')
        cubit.cmd('del vol all')
        cubit.cmd('set info on')
        cubit.cmd('set echo on')
        cubit.cmd('set journal on')
        return

    if number_of_chunks > 1:
        factor, minvalue, maxvalue, inv_length = prepare_equivalence_new()
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        if starting_tolerance:
            tol = starting_tolerance
        else:
            tol = minvalue / 20.
        if step_tolerance:
            step_tol = step_tolerance
        else:
            step_tol = minvalue / 20.

        cubit.cmd('group \'coincident_lateral_nodes\' add \
                   Node in face in group lateral')

        isempty = False
        while not isempty:
            isempty = merging_node_new(tol, clean=True, graphic_debug=False)
            tol = tol + step_tol
            if tol > maxvalue * 1.5:
                raise MergingError(
                    'tolerance greater than the max length of the edges, \
                    please check the mesh')

    # if checknodes and checklines:
    for ig, g in enumerate(store_group_name):
        cubit.cmd('block ' + str(ig + 1) + ' hex in group ' + g)
        cubit.cmd('block ' + str(ig + 1) + ' name "vol' + str(ig + 1) + '"')
        print('block ' + str(ig + 1) + ' hex in group ' + g)
    for ig, g in enumerate(side_name):
        cubit.cmd('block ' + side_val[ig] + ' face in group ' + g)
        print('block ' + side_val[ig] + ' face in group ' + g)
        cubit.cmd('block ' + side_val[ig] +
                  ' name "' + side_block_name[ig] + '"')
    cubit.cmd('del group all')

    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    cubit.cmd('set journal on')
    # cubit.cmd('set error on')

##########################################################################################
#
# deprecated methods
#
##########################################################################################

def add_sea_layer(block=1001, optionsea=False):
    if optionsea:
            # sea=optionsea['sea']
        seaup = optionsea['seaup']
        sealevel = optionsea['sealevel']
        seathres = optionsea['seathres']
    else:
        # sea=False
        seaup = False
        sealevel = False
        seathres = False

    # TODO
    # add sea hex
    # change hex absoorbing....
    block_list = cubit.get_block_id_list()
    id_block = max(block for block in block_list if block < 1000)
    cubit.cmd('delete block ' + str(id_block))
    # sea
    command = 'block ' + str(id_block) + ' hex in node in face in block ' + \
        str(block) + ' with Z_coord < ' + str(seathres)
    cubit.cmd(command)
    command = "block " + str(id_block) + " name 'sea'"
    cubit.cmd(command)
    if not seaup:
        id_block += 1
        cmd = 'block ' + str(id_block) + ' hex in node in face in block ' +\
              str(block) + ' with (Z_coord > ' + str(seathres) + \
              ' and Z_coord < ' + str(sealevel) + ')'
        cubit.cmd(cmd)
        command = "block " + str(id_block) + " name 'shwater'"
        cubit.cmd(command)
    id_block += 1
    command = 'block ' + str(id_block) + ' hex in node in face in block ' + \
        str(block) + ' with Z_coord >= ' + str(sealevel)
    cubit.cmd(command)
    command = "block " + str(id_block) + " name 'continent'"
    cubit.cmd(command)


def importing_cubfiles_old(cubfiles):
    import re
    rule_st = re.compile("(.+)_[0-9]+\.")
    rule_ex = re.compile(".+_[0-9]+\.(.+)")
    rule_int = re.compile(".+_([0-9]+)\.")
    filenames = glob.glob(cubfiles)
    try:
        st = rule_st.findall(filenames[0])[0]
        ex = rule_ex.findall(filenames[0])[0]
        listflag = True
    except:
        ex = ''
        listflag = False
    if ex == 'cub':
        cubflag = True
    else:
        cubflag = False
    list_int = []
    fs = []
    try:
        for f in filenames:
            i = int(rule_int.findall(f)[0])
            list_int.append(i)
        list_int.sort()
        for i, ind in enumerate(list_int):
            f = st + '_' + str(ind) + '.' + ex
            fs.append(f)
    except:
        pass
    if listflag:
        filenames = fs
    else:
        pass
    return len(filenames), list_int, filenames, cubflag


def refine_closecurve(block=1001, closed_filenames=None, acis=True):
    from utilities import load_curves
    from boundary_definition import build_block_side, define_surf
    from mesh_volume import refine_inside_curve
    #
    #
    curves = []
    if not isinstance(closed_filenames, list):
        closed_filenames = [closed_filenames]
    for f in closed_filenames:
        print(f)
        if acis:
            curves = curves + load_curves(f)
    print(curves)
    blist = list(cubit.get_block_id_list())
    try:
        blist.remove(1001)
    except:
        pass
    try:
        blist.remove(1002)
    except:
        pass
    try:
        blist.remove(1003)
    except:
        pass
    try:
        blist.remove(1004)
    except:
        pass
    try:
        blist.remove(1005)
    except:
        pass
    try:
        blist.remove(1006)
    except:
        pass
    id_top = max(blist)
    cmd = 'group "coi" add node in hex in block ' + str(id_top)
    cubit.cmd(cmd)
    #
    id_inside_arc = None
    for c in map(int, curves[0].split()):   # curves is a list of one string
        c1001 = cubit.get_exodus_element_count(1001, "block")
        c1002 = cubit.get_exodus_element_count(1002, "block")
        c1003 = cubit.get_exodus_element_count(1003, "block")
        c1004 = cubit.get_exodus_element_count(1004, "block")
        c1005 = cubit.get_exodus_element_count(1005, "block")
        c1006 = cubit.get_exodus_element_count(1006, "block")
        #
        refine_inside_curve(c, ntimes=1, depth=1, block=block, surface=False)
        blist = list(cubit.get_block_id_list())
        cmd = 'create mesh geometry hex all except hex in \
              block all feature_angle 135'
        cubit.cmd(cmd)
        blist_after = list(cubit.get_block_id_list())
        [blist_after.remove(x) for x in blist]
        id_inside = max(blist_after)
        cmd = 'group "coi" add node in hex in block ' + str(id_inside)
        cubit.cmd(cmd)
        if id_inside_arc:
            cmd = 'del block ' + str(id_inside - 1)
            cubit.cmd(cmd)
        cmd = 'block ' + str(id_inside) + ' name "refined"'
        cubit.cmd(cmd)
        id_inside_arc = id_inside
        #
        _, _, _, _, _, top_surf, bottom_surf, surf_xmin, \
            surf_ymin, surf_xmax, surf_ymax = define_surf()
        #
        c1001_after = cubit.get_exodus_element_count(1001, "block")
        c1002_after = cubit.get_exodus_element_count(1002, "block")
        c1003_after = cubit.get_exodus_element_count(1003, "block")
        c1004_after = cubit.get_exodus_element_count(1004, "block")
        c1005_after = cubit.get_exodus_element_count(1005, "block")
        c1006_after = cubit.get_exodus_element_count(1006, "block")
        entity = 'face'
        if c1001_after != c1001:
            refname = entity + '_topo'
            build_block_side(top_surf, refname, obj=entity, id_0=1001)
        #
        if c1002_after != c1002:
            refname = entity + '_bottom'
            build_block_side(bottom_surf, refname, obj=entity, id_0=1002)
        #
        if c1003_after != c1003:
            refname = entity + '_abs_xmin'
            build_block_side(surf_xmin, refname, obj=entity, id_0=1003)
        #
        if c1004_after != c1004:
            refname = entity + '_abs_ymin'
            build_block_side(surf_ymin, refname, obj=entity, id_0=1004)
        #
        if c1005_after != c1005:
            refname = entity + '_abs_xmax'
            build_block_side(surf_xmax, refname, obj=entity, id_0=1005)
        #
        if c1006_after != c1006:
            refname = entity + '_abs_ymax'
            build_block_side(surf_ymax, refname, obj=entity, id_0=1006)
        #
        cmd = 'disassociate mesh from volume all'
        cubit.cmd(cmd)
        cmd = 'group "coi" add node in face in \
              block 1001 1002 1003 1004 1005 1006'
        cubit.cmd(cmd)
        cubit.cmd('del vol all')
        cubit.cmd('group "removedouble" add hex all except hex in block all')
        cubit.cmd('delete hex in removedouble')
        cubit.cmd('delet group removedouble')
        cmd = 'equivalence node in group coi tolerance 20'
        cubit.cmd(cmd)
    cmd = 'equivalence node all tolerance 10'
    cubit.cmd(cmd)
    cubit.cmd('del curve ' + ' '.join(str(x) for x in curves))


def collecting_merging(cpuxmin=0, cpuxmax=1, cpuymin=0, cpuymax=1, cpux=1,
                       cpuy=1, cubfiles=False, ckbound_method1=False,
                       ckbound_method2=False, merge_tolerance=None,
                       decimate=False):

    boundary_dict = {}
    ##
    try:
        from boundary_definition import check_bc, map_boundary
    except:
        pass
    #
    xmin, xmax, ymin, ymax, listfull = map_boundary(
        cpuxmin, cpuxmax, cpuymin, cpuymax, cpux, cpuy)
    #
    if cubfiles:
        nf, listiproc, filenames, cubflag = importing_cubfiles(cubfiles)
    else:
        nf = 0
        filenames = []
        iproc = 0
    #
    if nf > 0:
        for iproc, filename in zip(listiproc, filenames):
            try:
                if iproc in listfull:
                    if cubflag:
                        cubit.cmd('import cubit "' + filename + '"')
                    else:
                        cubit.cmd('import mesh geometry "' + filename +
                                  '" block all use nodeset sideset \
                                  feature_angle 135.00 linear merge')
                    if decimate:
                        cubit.cmd(
                            'refine volume all numsplit 1 bias 1.0 depth 1 ')
                    boundary = check_bc(iproc, xmin, xmax, ymin, ymax,
                                        cpux, cpuy, cpuxmin, cpuxmax,
                                        cpuymin, cpuymax)
                    boundary_dict[iproc] = boundary
                    list_vol = list(cubit.parse_cubit_list('volume', 'all'))
                    for v in list_vol:
                        cubit.cmd("disassociate mesh from volume " + str(v))
                        command = "del vol " + str(v)
                        cubit.cmd(command)
            except:
                cubit.cmd('import mesh geometry "' + filename +
                          '" block all use nodeset sideset \
                          feature_angle 135.00 linear merge')
                if decimate:
                    cubit.cmd('refine volume all numsplit 1 bias 1.0 depth 1 ')
                iproc = 0
                boundary = check_bc(iproc, xmin, xmax, ymin, ymax,
                                    cpux, cpuy, cpuxmin, cpuxmax,
                                    cpuymin, cpuymax)
                boundary_dict[iproc] = boundary
                list_vol = list(cubit.parse_cubit_list('volume', 'all'))
                for v in list_vol:
                    cubit.cmd("disassociate mesh from volume " + str(v))
                    command = "del vol " + str(v)
                    cubit.cmd(command)
        cubit.cmd('export mesh "tmp_collect_NOmerging.e" \
                  dimension 3 block all overwrite')
    else:
        if decimate:
            cubit.cmd('refine volume all numsplit 1 bias 1.0 depth 1 ')
        boundary = check_bc(iproc, xmin, xmax, ymin, ymax, cpux,
                            cpuy, cpuxmin, cpuxmax, cpuymin, cpuymax)
    #
    #
    # print(boundary_dict)
    block_list = cubit.get_block_id_list()
    for block in block_list:
        ty = cubit.get_block_element_type(block)
        if ty == 'HEX8':
            cubit.cmd('block ' + str(block) + ' name "vol' + str(block) + '"')
    #
    #
    print('chbound', ckbound_method1, ckbound_method2)

    if ckbound_method1 and not ckbound_method2 and len(filenames) != 1:
        # use the equivalence method for groups
        if isinstance(merge_tolerance, list):
            tol = merge_tolerance[0]
        elif merge_tolerance:
            tol = merge_tolerance
        else:
            tol = 100000
        #
        idiag = None
        # cubit.cmd('set info off')
        # cubit.cmd('set journal off')
        # cubit.cmd('set echo off')
        ind = 0
        for ix in range(cpuxmin, cpuxmax):
            for iy in range(cpuymin, cpuymax):
                ind = ind + 1
                iproc = iy * cpux + ix
                print('******************* ', iproc, ind, '/', len(listfull))
                #
                #   ileft    |   iproc
                #  --------------------
                #   idiag    |   idown
                #
                #
                if iproc not in xmin and iproc not in ymin:
                    ileft = iy * cpux + ix - 1
                    idown = (iy - 1) * cpux + ix
                    idiag = idown - 1
                elif iproc in xmin and iproc in ymin:
                    ileft = iproc
                    idown = iproc
                    idiag = None
                elif iproc in xmin:
                    ileft = iproc
                    idown = (iy - 1) * cpux + ix
                    idiag = idown
                elif iproc in ymin:
                    ileft = iy * cpux + ix - 1
                    idown = iproc
                    idiag = ileft
                #
                print(iproc, ileft, idiag, idown)
                if iproc != idown:
                    nup = boundary_dict[iproc]['nodes_surf_ymin']
                    ndow = boundary_dict[idown]['nodes_surf_ymax']
                    merge_node_ck(nup, ndow)

                    if idiag != idown:
                        if iproc in ymax and iproc not in xmin:
                            # node in curve chunck left up... r u
                            nlu = boundary_dict[iproc]['node_curve_xminymax']
                            nru = boundary_dict[ileft]['node_curve_xmaxymax']
                            merge_node(nlu, nru)
                        if iproc in xmax:
                            # node in curve chunck left up... r u
                            nrd = boundary_dict[iproc]['node_curve_xmaxymin']
                            nru = boundary_dict[idown]['node_curve_xmaxymax']
                            merge_node(nrd, nru)
                        # node in curve chunck right up... r u
                        nru = boundary_dict[iproc]['node_curve_xminymin']
                        nrd = boundary_dict[idown]['node_curve_xminymax']
                        nld = boundary_dict[idiag]['node_curve_xmaxymax']
                        nlu = boundary_dict[ileft]['node_curve_xmaxymin']
                        merge_node_4(nru, nrd, nld, nlu)
                    elif iproc in xmin:
                        # node in curve chunck right up... r u
                        nlu = boundary_dict[iproc]['node_curve_xminymin']
                        nld = boundary_dict[idown]['node_curve_xminymax']
                        merge_node(nld, nlu)
                        # node in curve chunck right up... r u
                        nru = boundary_dict[iproc]['node_curve_xmaxymin']
                        nrd = boundary_dict[idown]['node_curve_xmaxymax']
                        merge_node(nrd, nru)

                #
                if iproc != ileft:
                    nright = boundary_dict[iproc]['nodes_surf_xmin']
                    nleft = boundary_dict[ileft]['nodes_surf_xmax']
                    merge_node_ck(nright, nleft)
                    #
                    #
                    if iproc in ymin:
                        # node in curve chunck right down... r u
                        nrd = boundary_dict[iproc]['node_curve_xminymin']
                        nld = boundary_dict[ileft]['node_curve_xmaxymin']
                        merge_node(nrd, nld)
                    if iproc in ymax:
                        # node in curve chunck right up... r u
                        nru = boundary_dict[iproc]['node_curve_xminymax']
                        nlu = boundary_dict[ileft]['node_curve_xmaxymax']
                        merge_node(nlu, nru)

        cubit.cmd('set info on')
        cubit.cmd('set echo on')
        cubit.cmd('set journal on')

        #
        #
        cmd = 'group "negativejac" add quality hex all Jacobian high'
        cubit.cmd(cmd)
        group_id_1 = cubit.get_id_from_name("negativejac")
        n1 = cubit.get_group_nodes(group_id_1)
        if len(n1) != 0:
            print('error, negative jacobian after the equivalence node command, \
                  use --merge2 instead of --equivalence/--merge/--merge1')
    elif ckbound_method2 and not ckbound_method1 and len(filenames) != 1:
        if isinstance(merge_tolerance, list):
            tol = merge_tolerance[0]
        elif merge_tolerance:
            tol = merge_tolerance
        else:
            tol = 100000
        #
        idiag = None
        for ix in range(cpuxmin, cpuxmax):
            for iy in range(cpuymin, cpuymax):
                iproc = iy * cpux + ix
                print('******************* ', iproc)
                #
                #   ileft    |   iproc
                #  --------------------
                #   idiag    |   idown
                #
                #
                if iproc not in xmin and iproc not in ymin:
                    ileft = iy * cpux + ix - 1
                    idown = (iy - 1) * cpux + ix
                    idiag = idown - 1
                elif iproc in xmin and iproc in ymin:
                    ileft = iproc
                    idown = iproc
                elif iproc in xmin:
                    ileft = iproc
                    idown = (iy - 1) * cpux + ix
                    idiag = idown
                elif iproc in ymin:
                    ileft = iy * cpux + ix - 1
                    idown = iproc
                    idiag = ileft
                #
                #
                if iproc != idown:
                    nup = boundary_dict[iproc]['nodes_surf_ymin']
                    ndow = boundary_dict[idown]['nodes_surf_ymax']
                    for n1, n2 in zip(nup, ndow):
                        cubit.cmd('equivalence node ' + str(n1) +
                                  ' ' + str(n2) + ' tolerance ' + str(tol))
                    if idiag != idown:
                        if iproc in ymax and iproc not in xmin:
                            # node in curve chunck left up... r u
                            nlu = boundary_dict[iproc]['node_curve_xminymax']
                            nru = boundary_dict[ileft]['node_curve_xmaxymax']
                            for n in zip(nlu, nru):
                                cubit.cmd('equivalence node ' +
                                          ' '.join(str(x) for x in n) +
                                          ' tolerance ' + str(tol))

                        # node in curve chunck right up... r u
                        nru = boundary_dict[iproc]['node_curve_xminymin']
                        nrd = boundary_dict[idown]['node_curve_xminymax']
                        nld = boundary_dict[idiag]['node_curve_xmaxymax']
                        nlu = boundary_dict[ileft]['node_curve_xmaxymin']
                        for n in zip(nru, nrd, nlu, nld):
                            cubit.cmd('equivalence node ' +
                                      ' '.join(str(x) for x in n) +
                                      ' tolerance ' + str(tol))

                    elif iproc in xmin:
                        # node in curve chunck right up... r u
                        nru = boundary_dict[iproc]['node_curve_xminymin']
                        nrd = boundary_dict[idown]['node_curve_xminymax']
                        for n in zip(nru, nrd):
                            cubit.cmd('equivalence node ' +
                                      ' '.join(str(x) for x in n) +
                                      ' tolerance ' + str(tol))
                #
                #
                if iproc != ileft:
                    nright = boundary_dict[iproc]['nodes_surf_xmin']
                    nleft = boundary_dict[ileft]['nodes_surf_xmax']
                    for n1, n2 in zip(nleft, nright):
                        cubit.cmd('equivalence node ' + str(n1) +
                                  ' ' + str(n2) + ' tolerance ' + str(tol))
                    #
                    #
                    if iproc in ymin:
                        # node in curve chunck right down... r u
                        nrd = boundary_dict[iproc]['node_curve_xminymin']
                        nld = boundary_dict[ileft]['node_curve_xmaxymin']
                        for n in zip(nrd, nld):
                            cubit.cmd('equivalence node ' +
                                      ' '.join(str(x) for x in n) +
                                      ' tolerance ' + str(tol))
                    if iproc in ymax:
                        # node in curve chunck right up... r u
                        nru = boundary_dict[iproc]['node_curve_xminymax']
                        nlu = boundary_dict[ileft]['node_curve_xmaxymax']
                        for n in zip(nru, nlu):
                            cubit.cmd('equivalence node ' +
                                      ' '.join(str(x) for x in n) +
                                      ' tolerance ' + str(tol))
        #
        #
        cmd = 'topology check coincident node face all tolerance ' + \
            str(tol * 2) + ' nodraw brief result group "checkcoinc"'
        cubit.silent_cmd(cmd)
        group_id_1 = cubit.get_id_from_name("checkcoinc")
        if group_id_1 != 0:
            n1 = cubit.get_group_nodes(group_id_1)
            if len(n1) != 0:
                print('error, coincident nodes after the equivalence \
                      node command, check the tolerance')
                import sys
                sys.exit()
        cmd = 'group "negativejac" add quality hex all Jacobian high'
        cubit.cmd(cmd)
        group_id_1 = cubit.get_id_from_name("negativejac")
        n1 = cubit.get_group_nodes(group_id_1)
        if len(n1) != 0:
            print('error, negative jacobian after the equivalence node command, \
                  check the mesh')
    elif ckbound_method1 and ckbound_method2 and len(filenames) != 1:
        block_list = cubit.get_block_id_list()
        i = -1
        for block in block_list:
            ty = cubit.get_block_element_type(block)
            if ty == 'HEX8':
                i = i + 1
                if isinstance(merge_tolerance, list):
                    try:
                        tol = merge_tolerance[i]
                    except:
                        tol = merge_tolerance[-1]
                elif merge_tolerance:
                    tol = merge_tolerance
                else:
                    tol = 1
                cmd = 'topology check coincident node face in hex in block ' +\
                    str(block) + ' tolerance ' + str(tol) + \
                    ' nodraw brief result group "b' + str(block) + '"'
                cubit.cmd(cmd)
                print(cmd)
                cmd = 'equivalence node in group b' + \
                    str(block) + ' tolerance ' + str(tol)
                cubit.cmd(cmd)
                print(cmd)
        if isinstance(merge_tolerance, list):
            tol = max(merge_tolerance)
        elif merge_tolerance:
            tol = merge_tolerance
        else:
            tol = 1
        #
        #
        cmd = 'topology check coincident node face all tolerance ' +\
            str(tol) + ' nodraw brief result group "checkcoinc"'
        cubit.silent_cmd(cmd)
        group_id_1 = cubit.get_id_from_name("checkcoinc")
        if group_id_1 != 0:
            n1 = cubit.get_group_nodes(group_id_1)
            if len(n1) != 0:
                print('error, coincident nodes after the equivalence node \
                       command, check the tolerance')
                import sys
                sys.exit()
        cmd = 'group "negativejac" add quality hex all Jacobian high'
        cubit.silent_cmd(cmd)
        group_id_1 = cubit.get_id_from_name("negativejac")
        n1 = cubit.get_group_nodes(group_id_1)
        if len(n1) != 0:
            print('error, negative jacobian after the equivalence node command, \
                  use --merge instead of --equivalence')


def collect_old(cpuxmin=0, cpuxmax=1, cpuymin=0, cpuymax=1, cpux=1, cpuy=1,
                cubfiles=False, ckbound_method1=False, ckbound_method2=False,
                merge_tolerance=None, curverefining=False,
                outfilename='totalmesh_merged', qlog=False,
                export2SPECFEM3D=False, listblock=None,
                listflag=None, outdir='.', add_sea=False, decimate=False,
                cpml=False, cpml_size=False, top_absorbing=False, hex27=False):
    #
    # cubit.cmd('set journal error off')
    # cubit.cmd('set verbose error off')
    collecting_merging(cpuxmin, cpuxmax, cpuymin, cpuymax, cpux, cpuy,
                       cubfiles=cubfiles, ckbound_method1=ckbound_method1,
                       ckbound_method2=ckbound_method2,
                       merge_tolerance=merge_tolerance, decimate=decimate)
    # cubit.cmd('set journal error on')
    # cubit.cmd('set verbose error on')
    #
    if curverefining:
        block = 1001  # topography
        refine_closecurve(block, curverefining, acis=True)
    #
    #
    if add_sea:
        block = 1001
        add_sea_layer(block=block)

    outdir2 = '/'.join(x for x in outfilename.split('/')[:-1])
    if outdir2 == '':
        outdir2 = outdir + '/'
    else:
        outdir2 = outdir + '/' + outdir2 + '/'

    import os
    try:
        os.makedirs(outdir2)
    except OSError:
        pass

    cubit.cmd('compress all')
    command = "export mesh '" + outdir2 + outfilename + \
        ".e' block all overwrite xml '" + outdir2 + outfilename + ".xml'"
    cubit.cmd(command)
    f = open(outdir2 + 'blocks.dat', 'w')
    blocks = cubit.get_block_id_list()
    #
    for block in blocks:
        name = cubit.get_exodus_entity_name('block', block)
        element_count = cubit.get_exodus_element_count(block, "block")
        nattrib = cubit.get_block_attribute_count(block)
        attr = [cubit.get_block_attribute_value(
            block, x) for x in range(0, nattrib)]
        ty = cubit.get_block_element_type(block)
        f.write(str(block) + ' ; ' + name + ' ; nattr ' + str(nattrib) +
                ' ; ' + ' '.join(str(x) for x in attr) + ' ; ' + ty + ' ' +
                str(element_count) + '\n')
    f.close()
    #
    #
    cubit.cmd('set info echo journ off')
    cmd = 'del group all'
    cubit.silent_cmd(cmd)
    cubit.cmd('set info echo journ on')
    #
    command = "save as '" + outdir2 + outfilename + ".cub' overwrite"
    cubit.cmd(command)
    #
    print('end meshing')
    #
    #
    if qlog:
        print('\n\nQUALITY CHECK.... ***************\n\n')
        import quality_log
        tq = open(outdir2 + outfilename + '.quality', 'w')
        max_skewness, min_length = quality_log.quality_log(tq)
    #
    #
    #
    if export2SPECFEM3D:
        e2SEM(files=False, listblock=listblock,
              listflag=listflag, outdir=outdir,
              cpml=cpml, cpml_size=cpml_size,
              top_absorbing=top_absorbing, hex27=hex27)


def e2SEM_old(files=False, listblock=None, listflag=None, outdir='.',
              cpml=False, cpml_size=False, top_absorbing=False, hex27=False):
    import glob
    if files:
        filenames = glob.glob(files)
        for f in filenames:
            print(f)
            extension = f.split('.')[-1]
            if extension == 'cub':
                cubit.cmd('open "' + f + '"')
            elif extension == 'e':
                cubit.cmd('import mesh "' + f + '" no_geom')
            else:
                print(extension)
    if listblock and listflag:
        pass
    else:
        listblock = []
        listflag = []
        block_list = list(cubit.get_block_id_list())
        for block in block_list:
            ty = cubit.get_block_element_type(block)
            if 'HEX' in ty:
                listblock.append(block)
                # listflag.append(block)
        listflag = range(1, len(block_list) + 1)
    #
    for ib, iflag in zip(listblock, listflag):
        cubit.cmd("block " + str(ib) + " attribute count 1")
        cubit.cmd("block " + str(ib) + " attribute index 1 " + str(iflag))
    #
    import cubit2specfem3d
    import os
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    cubit2specfem3d.export2SPECFEM3D(outdir, cpml=cpml, cpml_size=cpml_size,
                                     top_absorbing=top_absorbing, hex27=hex27)


def prepare_equivalence(nodes1, nodes2):
    cubit.cmd('set info off')
    cubit.cmd('set echo off')
    cubit.cmd('set journal off')
    length = {}
    for ns in zip(nodes1, nodes2):
        cmd = 'group "tmpn" add edge in node ' + ' '.join(str(n) for n in ns)
        cubit.cmd(cmd)
        ge = cubit.get_id_from_name("tmpn")
        e1 = cubit.get_group_edges(ge)
        lengthmin = 1e9
        for e in e1:
            lengthmin = min(lengthmin, cubit.get_mesh_edge_length(e))
        length[ns] = lengthmin * .5
        cubit.cmd('delete group ' + str(ge))
    minvalue = min(length.values())
    maxvalue = max(length.values())
    print('min lentgh: ', minvalue, 'max lentgh: ', maxvalue)
    nbin = int((maxvalue / minvalue) / 2.) + 1
    factor = (maxvalue - minvalue) / nbin
    dic_new = {}
    for k in length.keys():
        if factor != 0.:
            dic_new[k] = int((length[k] - minvalue) / factor)
        else:
            dic_new[k] = 0.
    inv_length = invert_dict(dic_new)
    print(inv_length.keys(), factor, minvalue)
    ks = inv_length.keys()
    ks.sort()
    for k in range(0, len(inv_length.keys()) - 1):
        inv_length[ks[k]] = inv_length[ks[k]] + inv_length[ks[k + 1]]
    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    cubit.cmd('set journal on')
    return factor, minvalue, inv_length


def merge_node_ck(n1, n2):
    factor, minvalue, inv_length = prepare_equivalence(n1, n2)

    cubit.cmd('set info off')
    cubit.cmd('set echo off')
    cubit.cmd('set journal off')
    # cubit.cmd('set error off')

    for k in inv_length.keys()[:-1]:
        if len(inv_length[k]) > 0:
            cmd = 'equivalence node ' + \
                  ' '.join(' '.join(str(n) for n in x)
                           for x in inv_length[k]) +\
                  ' tolerance ' + str(k * factor + minvalue / 3.)
            cubit.cmd(cmd)
            print('equivalence ' + str(len(inv_length[k])) +\
                  ' couples of nodes -  tolerance ' + \
                  str(k * factor + minvalue / 3.))

    cubit.cmd('group "checkmerge" add node ' +
              ' '.join(str(n) for n in n1) +
              ' ' + ' '.join(str(n) for n in n2))

    idg = cubit.get_id_from_name('checkmerge')
    remainnodes = cubit.get_group_nodes(idg)
    print('from ' + str(len(n1) + len(n2)) + ' nodes -> ' + \
          str(len(remainnodes)) + ' nodes')
    if len(n1) != len(remainnodes):
        print('equivalence ' + str(len(remainnodes)) + \
              ' couples of nodes -  tolerance ' + str(minvalue / 3.))
        cubit.cmd('set info on')
        cubit.cmd('set echo on')
        cubit.cmd('set journal on')
        cmd = 'equivalence node in group ' + \
            str(idg) + ' tolerance ' + str(minvalue / 3.)
        cubit.cmd(cmd)
        cmd = 'block 3000 node in group ' + str(idg)
        cubit.cmd(cmd)

    if len(n1) != len(remainnodes):
        cubit.cmd('export mesh "error_merging.e" \
                  dimension 3 block all overwrite')
        cubit.cmd('save as "error_merging.cub" \
                  dimension 3 block all overwrite')
        print('error merging ')
        if False:
            import sys
            sys.exit(2)

    cubit.cmd('delete group checkmerge')
    cubit.cmd('delete block 3000')

    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    cubit.cmd('set journal on')


def merge_node(n1, n2):
    factor, minvalue, inv_length = prepare_equivalence(n1, n2)

    cubit.cmd('set info off')
    cubit.cmd('set echo off')
    cubit.cmd('set journal off')

    for k in inv_length.keys()[:-1]:
        if len(inv_length[k]) > 0:
            cmd = 'equivalence node ' + \
                  ' '.join(' '.join(str(n) for n in x)
                           for x in inv_length[k]) +\
                  ' tolerance ' + str(k * factor + minvalue / 3.)

            cubit.cmd(cmd)
            print('equivalence ' + str(len(inv_length[k])) + \
                  ' couples of nodes -  tolerance ' + \
                  str(k * factor + minvalue / 3.))

    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    cubit.cmd('set journal on')


def prepare_equivalence_4(nodes1, nodes2, nodes3, nodes4):
    cubit.cmd('set info off')
    cubit.cmd('set echo off')
    cubit.cmd('set journal off')
    length = {}
    nodes = [nodes1, nodes2, nodes3, nodes4]
    check = map(len, nodes)
    checked_nodes = []
    for ind, iflag in enumerate(check):
        if iflag:
            checked_nodes = checked_nodes + nodes[ind]

    cmd = 'group "tmpn" add edge in node ' + \
        ' '.join(str(n) for n in checked_nodes)
    cubit.cmd(cmd)
    ge = cubit.get_id_from_name("tmpn")
    e1 = cubit.get_group_edges(ge)
    lengthmin = 1e9
    for e in e1:
        lengthmin = min(lengthmin, cubit.get_mesh_edge_length(e))
        length[e] = lengthmin * .5
    cubit.cmd('delete group ' + str(ge))
    try:
        minvalue = min(length.values())
        maxvalue = max(length.values())
    except:
        try:
            print(nodes)
            print('edges ', e1)
        except:
            pass
        minvalue = 10.
        maxvalue = 2000.
    print('min lentgh: ', minvalue, 'max lentgh: ', maxvalue)
    nbin = int((maxvalue / minvalue) / 2.) + 1
    factor = (maxvalue - minvalue) / nbin
    dic_new = {}
    for k in length.keys():
        if factor != 0.:
            dic_new[k] = int((length[k] - minvalue) / factor)
        else:
            dic_new[k] = 0.
    inv_length = invert_dict(dic_new)
    print(inv_length.keys(), factor, minvalue)
    ks = inv_length.keys()
    ks.sort()
    for k in range(0, len(inv_length.keys()) - 1):
        inv_length[ks[k]] = inv_length[ks[k]] + inv_length[ks[k + 1]]
    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    cubit.cmd('set journal on')
    return factor, minvalue, inv_length


def ording_z(nodes):
    def get_z(node):
        x, y, z = cubit.get_nodal_coordinates(node)
        return z
    d = [(get_z(node), node) for node in nodes]
    d.sort()
    return [x[1] for x in d]


def merge_node_4(n1, n2, n3, n4, newmethod=True):
    if newmethod:
        print("merge node 4 side")
        n1o = ording_z(n1)
        n2o = ording_z(n2)
        n3o = ording_z(n3)
        n4o = ording_z(n4)
        for ln in zip(n1o, n2o, n3o, n4o):
            cmd = 'equivalence node ' + \
                ' '.join(str(n) for n in ln) + ' tolerance 10000 '
            cubit.cmd(cmd)

    else:
        factor, minvalue, inv_length = prepare_equivalence_4(n1, n2, n3, n4)

        for k in inv_length.keys()[:-1]:
            if len(inv_length[k]) > 1:
                try:
                    for x in inv_length[k]:
                        if type(x) is not list:
                            x = [x]
                        else:
                            pass
                    cmd = 'equivalence node ' + \
                        ' '.join(' '.join(str(n) for n in x)) + \
                        ' tolerance ' + str(k * factor + minvalue / 3.)
                except:
                    print(k, "***************************************** s")
                    print(inv_length[k])

                cubit.cmd(cmd)
                print('equivalence ' + str(len(inv_length[k])) +\
                      ' couples of nodes -  tolerance ' + \
                      str(k * factor + minvalue / 3.))
            if len(inv_length[k]) == 1:
                cmd = 'equivalence node ' + \
                    ' '.join(' '.join(str(n) for n in inv_length[k])) + \
                    ' tolerance ' + str(k * factor + minvalue / 3.)
                cubit.cmd(cmd)
                print('equivalence ' + str(len(inv_length[k])) + \
                      ' couples of nodes -  tolerance ' + \
                      str(k * factor + minvalue / 3.))

##########################################################################################
#
# end deprecated methods
#
##########################################################################################

collect = collect_new
define_blocks = collect_new

