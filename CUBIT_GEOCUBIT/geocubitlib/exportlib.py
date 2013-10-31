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
# GEOCUBIT is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation, either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# GEOCUBIT is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with GEOCUBIT.  If not, see <http://www.gnu.org/licenses/>.         #
#                                                                           #
#############################################################################
#
try:
    import start as start
    cubit                   = start.start_cubit()
except:
    try:
        import cubit
    except:
        print 'error importing cubit, check if cubit is installed'
        pass

import glob

def add_sea_layer(block=1001,optionsea=False):
    if optionsea:
            sea=optionsea['sea']
            seaup=optionsea['seaup']
            sealevel=optionsea['sealevel']
            seathres=optionsea['seathres']
    else:
            sea=False
            seaup=False
            sealevel=False
            seathres=False
    
    ######TODO
    #add sea hex
    #change hex absoorbing....
    block_list=cubit.get_block_id_list()
    id_block = max(block for block in block_list if block<1000)
    cubit.cmd('delete block '+str(id_block))
    #sea
    command= 'block '+str(id_block)+' hex in node in face in block '+str(block)+' with Z_coord < '+str(seathres)
    cubit.cmd(command)
    command = "block "+str(id_block)+" name 'sea'"
    cubit.cmd(command)
    if not seaup:
        id_block+=1
        command= 'block '+str(id_block)+' hex in node in face in block '+str(block)+' with (Z_coord > '+str(seathres)+' and Z_coord < '+str(sealevel)+')'
        cubit.cmd(command)
        command = "block "+str(id_block)+" name 'shwater'"
        cubit.cmd(command)
    id_block+=1
    command= 'block '+str(id_block)+' hex in node in face in block '+str(block)+' with Z_coord >= '+str(sealevel)
    cubit.cmd(command)
    command = "block "+str(id_block)+" name 'continent'"
    cubit.cmd(command)
    
    
def importing_cubfiles(cubfiles):
    import re
    rule_st=re.compile("(.+)_[0-9]+\.")
    rule_ex=re.compile(".+_[0-9]+\.(.+)")
    rule_int=re.compile(".+_([0-9]+)\.")
    filenames=glob.glob(cubfiles)
    try:
        st = rule_st.findall(filenames[0])[0]
        ex = rule_ex.findall(filenames[0])[0]
        listflag=True
    except:
        ex=''
        listflag=False
    if ex == 'cub':
        cubflag=True
    else:
        cubflag=False
    list_int=[]
    fs=[]
    try:
        for f in filenames:
            i=int(rule_int.findall(f)[0])
            list_int.append(i)
        list_int.sort()
        for i,ind in enumerate(list_int):
            f=st+'_'+str(ind)+'.'+ex
            fs.append(f)
    except:
        pass
    if listflag:
        filenames=fs
    else:
        pass
    return len(filenames),list_int,filenames,cubflag
    
    
    
    
    



def refine_closecurve(block=1001,closed_filenames=None,acis=True):
    from utilities import load_curves
    from boundary_definition import build_block_side,define_surf
    from mesh_volume import refine_inside_curve
    #
    #
    curves=[]
    if not isinstance(closed_filenames,list): closed_filenames=[closed_filenames]
    for f in closed_filenames:
        print f
        if acis: 
            curves=curves+load_curves(f)
    print curves
    blist=list(cubit.get_block_id_list())
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
    id_top=max(blist)
    cmd='group "coi" add node in hex in block '+str(id_top)
    cubit.cmd(cmd)
    #
    id_inside_arc=None
    for c in map(int,curves[0].split()):   #curves is a list of one string
        c1001 = cubit.get_exodus_element_count(1001, "block")
        c1002 = cubit.get_exodus_element_count(1002, "block")
        c1003 = cubit.get_exodus_element_count(1003, "block")
        c1004 = cubit.get_exodus_element_count(1004, "block")
        c1005 = cubit.get_exodus_element_count(1005, "block")
        c1006 = cubit.get_exodus_element_count(1006, "block")
        #
        refine_inside_curve(c,ntimes=1,depth=1,block=block,surface=False)
        blist=list(cubit.get_block_id_list())
        cmd='create mesh geometry hex all except hex in block all feature_angle 135'
        cubit.cmd(cmd)
        blist_after=list(cubit.get_block_id_list())
        [blist_after.remove(x) for x in blist]
        id_inside=max(blist_after)
        cmd='group "coi" add node in hex in block '+str(id_inside)
        cubit.cmd(cmd)
        if id_inside_arc: 
                cmd='del block '+str(id_inside-1)
                cubit.cmd(cmd)
        cmd='block '+str(id_inside)+' name "refined"'
        cubit.cmd(cmd)
        id_inside_arc=id_inside                                            
        #
        _,_,_,_,_,top_surf,bottom_surf,surf_xmin,surf_ymin,surf_xmax,surf_ymax=define_surf()
        #
        c1001_after = cubit.get_exodus_element_count(1001, "block")
        c1002_after = cubit.get_exodus_element_count(1002, "block")
        c1003_after = cubit.get_exodus_element_count(1003, "block")
        c1004_after = cubit.get_exodus_element_count(1004, "block")
        c1005_after = cubit.get_exodus_element_count(1005, "block")
        c1006_after = cubit.get_exodus_element_count(1006, "block")
        entity='face'
        if c1001_after != c1001:
            refname=entity+'_topo'
            build_block_side(top_surf,refname,obj=entity,id_0=1001)
        #
        if c1002_after != c1002:
            refname=entity+'_bottom'
            build_block_side(bottom_surf,refname,obj=entity,id_0=1002)
        #
        if c1003_after != c1003:
            refname=entity+'_abs_xmin'
            build_block_side(surf_xmin,refname,obj=entity,id_0=1003)
        #                                                              
        if c1004_after != c1004:
            refname=entity+'_abs_ymin'
            build_block_side(surf_ymin,refname,obj=entity,id_0=1004)
        #                                                              
        if c1005_after != c1005:  
            refname=entity+'_abs_xmax'                                     
            build_block_side(surf_xmax,refname,obj=entity,id_0=1005)
        #
        if c1006_after != c1006:    
            refname=entity+'_abs_ymax'                                   
            build_block_side(surf_ymax,refname,obj=entity,id_0=1006)
        #
        cmd='disassociate mesh from volume all'
        cubit.cmd(cmd)
        cmd='group "coi" add node in face in block 1001 1002 1003 1004 1005 1006'
        cubit.cmd(cmd)
        cubit.cmd('del vol all')
        cubit.cmd('group "removedouble" add hex all except hex in block all')
        cubit.cmd('delete hex in removedouble')
        cubit.cmd('delet group removedouble')
        cmd='equivalence node in group coi tolerance 20'
        cubit.cmd(cmd)
    cmd='equivalence node all tolerance 10'
    cubit.cmd(cmd)
    cubit.cmd('del curve '+' '.join(str(x) for x in curves) )




def collecting_merging(cpuxmin=0,cpuxmax=1,cpuymin=0,cpuymax=1,cpux=1,cpuy=1,cubfiles=False,ckbound_method1=False,ckbound_method2=False,merge_tolerance=None,decimate=False):
    import glob
    import re
    #
    rule_st=re.compile("(.+)_[0-9]+\.")
    rule_ex=re.compile(".+_[0-9]+\.(.+)")
    rule_int=re.compile(".+_([0-9]+)\.")
    boundary_dict={}
    ##
    try:
        from boundary_definition import check_bc, map_boundary
    except:
        pass
    #
    xmin,xmax,ymin,ymax,listfull=map_boundary(cpuxmin,cpuxmax,cpuymin,cpuymax,cpux,cpuy)
    #
    if cubfiles:
        nf,listip,filenames,cubflag=importing_cubfiles(cubfiles)
    else:
        nf=0
        filenames=[]
        ip=0
    #
    if nf > 0:
        for ip,filename in zip(listip,filenames):
            try:
                if ip in listfull:
                    if cubflag:
                        cubit.cmd('import cubit "'+filename+'"')
                    else:
                        cubit.cmd('import mesh geometry "'+filename+'" block all use nodeset sideset feature_angle 135.00 linear merge')
                    if decimate: cubit.cmd('refine volume all numsplit 1 bias 1.0 depth 1 ')
                    boundary=check_bc(ip,xmin,xmax,ymin,ymax,cpux,cpuy,cpuxmin,cpuxmax,cpuymin,cpuymax)
                    boundary_dict[ip]=boundary
                    list_vol=list(cubit.parse_cubit_list('volume','all'))
                    for v in list_vol:
                        cubit.cmd("disassociate mesh from volume "+str(v))
                        command = "del vol "+str(v)
                        cubit.cmd(command)
            except:
                cubit.cmd('import mesh geometry "'+filename+'" block all use nodeset sideset feature_angle 135.00 linear merge')
                if decimate: cubit.cmd('refine volume all numsplit 1 bias 1.0 depth 1 ')
                ip=0
                boundary=check_bc(ip,xmin,xmax,ymin,ymax,cpux,cpuy,cpuxmin,cpuxmax,cpuymin,cpuymax)
                boundary_dict[ip]=boundary
                list_vol=list(cubit.parse_cubit_list('volume','all'))
                for v in list_vol:
                    cubit.cmd("disassociate mesh from volume "+str(v))
                    command = "del vol "+str(v)
                    cubit.cmd(command)
        cubit.cmd('export mesh "tmp_collect_NOmerging.e" dimension 3 block all overwrite')
    else:
        if decimate: cubit.cmd('refine volume all numsplit 1 bias 1.0 depth 1 ')
        boundary=check_bc(ip,xmin,xmax,ymin,ymax,cpux,cpuy,cpuxmin,cpuxmax,cpuymin,cpuymax)
    #
    #
    #print boundary_dict
    block_list=cubit.get_block_id_list()
    for block in block_list:
        ty=cubit.get_block_element_type(block)
        if ty == 'HEX8':
            cubit.cmd('block '+str(block)+' name "vol'+str(block)+'"')
    #
    #
    print 'chbound',ckbound_method1,ckbound_method2
    
    
    if ckbound_method1 and not ckbound_method2 and len(filenames) != 1:
        #use the equivalence method for groups
        if isinstance(merge_tolerance,list):
            tol=merge_tolerance[0]
        elif merge_tolerance:
            tol=merge_tolerance
        else:
            tol=100000
        #
        idiag=None
        #cubit.cmd('set info off')
        #cubit.cmd('set journal off')
        #cubit.cmd('set echo off')
        ind=0
        for ix in range(cpuxmin,cpuxmax):
            for iy in range(cpuymin,cpuymax):
                ind=ind+1
                ip=iy*cpux+ix
                print '******************* ',ip, ind,'/',len(listfull)
                #
                #   ileft    |   ip
                #  --------------------
                #   idiag    |   idown
                #
                #
                if ip not in xmin and ip not in ymin:
                    ileft=iy*cpux+ix-1
                    idown=(iy-1)*cpux+ix
                    idiag=idown-1
                elif ip in xmin and ip in ymin:
                    ileft=ip
                    idown=ip
                    idiag=None
                elif ip in xmin:
                    ileft=ip
                    idown=(iy-1)*cpux+ix
                    idiag=idown
                elif ip in ymin:
                    ileft=iy*cpux+ix-1
                    idown=ip
                    idiag=ileft
                #
                print ip,ileft,idiag,idown
                if ip != idown:
                    nup=boundary_dict[ip]['nodes_surf_ymin']
                    ndow=boundary_dict[idown]['nodes_surf_ymax']
                    merge_node_ck(nup,ndow)
                 
                    if idiag != idown:
                        if ip in ymax and ip not in xmin:
                            nlu=boundary_dict[ip]['node_curve_xminymax'] #node in curve chunck left up... r u
                            nru=boundary_dict[ileft]['node_curve_xmaxymax']
                            merge_node(nlu,nru)
                        if ip in xmax:
                            nrd=boundary_dict[ip]['node_curve_xmaxymin'] #node in curve chunck left up... r u
                            nru=boundary_dict[idown]['node_curve_xmaxymax']
                            merge_node(nrd,nru)
                        nru=boundary_dict[ip]['node_curve_xminymin'] #node in curve chunck right up... r u
                        nrd=boundary_dict[idown]['node_curve_xminymax']
                        nld=boundary_dict[idiag]['node_curve_xmaxymax']
                        nlu=boundary_dict[ileft]['node_curve_xmaxymin']
                        merge_node_4(nru,nrd,nld,nlu)
                    elif ip in xmin:
                        nlu=boundary_dict[ip]['node_curve_xminymin'] #node in curve chunck right up... r u
                        nld=boundary_dict[idown]['node_curve_xminymax']
                        merge_node(nld,nlu)
                        nru=boundary_dict[ip]['node_curve_xmaxymin'] #node in curve chunck right up... r u
                        nrd=boundary_dict[idown]['node_curve_xmaxymax']
                        merge_node(nrd,nru)
                        
                        
                        
                        
                #
                if ip != ileft:
                    nright=boundary_dict[ip]['nodes_surf_xmin']
                    nleft=boundary_dict[ileft]['nodes_surf_xmax']
                    merge_node_ck(nright,nleft)
                    #
                    #
                    if ip in ymin:
                        nrd=boundary_dict[ip]['node_curve_xminymin'] #node in curve chunck right down... r u
                        nld=boundary_dict[ileft]['node_curve_xmaxymin']
                        merge_node(nrd,nld)
                    if ip in ymax:
                        nru=boundary_dict[ip]['node_curve_xminymax'] #node in curve chunck right up... r u
                        nlu=boundary_dict[ileft]['node_curve_xmaxymax']
                        merge_node(nlu,nru)
        
        cubit.cmd('set info on')
        cubit.cmd('set echo on')
        cubit.cmd('set journal on')
        
        
        #
        #
        cmd='group "negativejac" add quality hex all Jacobian high'
        cubit.cmd(cmd) 
        group_id_1=cubit.get_id_from_name("negativejac")
        n1=cubit.get_group_nodes(group_id_1)
        if len(n1) != 0:
            print 'error, negative jacobian after the equivalence node command, use --merge2 instead of --equivalence/--merge/--merge1'
    elif ckbound_method2 and not ckbound_method1 and len(filenames) != 1:
        if isinstance(merge_tolerance,list):
            tol=merge_tolerance[0]
        elif merge_tolerance:
            tol=merge_tolerance
        else:
            tol=100000
        #
        idiag=None
        for ix in range(cpuxmin,cpuxmax):
            for iy in range(cpuymin,cpuymax):
                ip=iy*cpux+ix
                print '******************* ',ip
                #
                #   ileft    |   ip
                #  --------------------
                #   idiag    |   idown
                #
                #
                if ip not in xmin and ip not in ymin:
                    ileft=iy*cpux+ix-1
                    idown=(iy-1)*cpux+ix
                    idiag=idown-1
                elif ip in xmin and ip in ymin:
                    ileft=ip
                    idown=ip
                elif ip in xmin:
                    ileft=ip
                    idown=(iy-1)*cpux+ix
                    idiag=idown
                elif ip in ymin:
                    ileft=iy*cpux+ix-1
                    idown=ip
                    idiag=ileft
                #
                #
                if ip != idown:
                    nup=boundary_dict[ip]['nodes_surf_ymin']
                    ndow=boundary_dict[idown]['nodes_surf_ymax']
                    for n1,n2 in zip(nup,ndow):
                        cubit.cmd('equivalence node '+str(n1)+' '+str(n2)+' tolerance '+str(tol))
                    if idiag != idown:
                        if ip in ymax and ip not in xmin:
                            nlu=boundary_dict[ip]['node_curve_xminymax'] #node in curve chunck left up... r u
                            nru=boundary_dict[ileft]['node_curve_xmaxymax']
                            for n in zip(nlu,nru):
                                cubit.cmd('equivalence node '+' '.join(str(x) for x in n)+' tolerance '+str(tol))
                        nru=boundary_dict[ip]['node_curve_xminymin'] #node in curve chunck right up... r u
                        nrd=boundary_dict[idown]['node_curve_xminymax']
                        nld=boundary_dict[idiag]['node_curve_xmaxymax']
                        nlu=boundary_dict[ileft]['node_curve_xmaxymin']
                        for n in zip(nru,nrd,nlu,nld):
                            cubit.cmd('equivalence node '+' '.join(str(x) for x in n)+' tolerance '+str(tol))
                    elif ip in xmin:
                        nru=boundary_dict[ip]['node_curve_xminymin'] #node in curve chunck right up... r u
                        nrd=boundary_dict[idown]['node_curve_xminymax']
                        for n in zip(nru,nrd):
                            cubit.cmd('equivalence node '+' '.join(str(x) for x in n)+' tolerance '+str(tol))
                #
                #
                if ip != ileft:
                    nright=boundary_dict[ip]['nodes_surf_xmin']
                    nleft=boundary_dict[ileft]['nodes_surf_xmax']
                    for n1,n2 in zip(nleft,nright):
                        cubit.cmd('equivalence node '+str(n1)+' '+str(n2)+' tolerance '+str(tol))
                    #
                    #
                    if ip in ymin:
                        nrd=boundary_dict[ip]['node_curve_xminymin'] #node in curve chunck right down... r u
                        nld=boundary_dict[ileft]['node_curve_xmaxymin']
                        for n in zip(nrd,nld):
                            cubit.cmd('equivalence node '+' '.join(str(x) for x in n)+' tolerance '+str(tol))
                    if ip in ymax:
                        nru=boundary_dict[ip]['node_curve_xminymax'] #node in curve chunck right up... r u
                        nlu=boundary_dict[ileft]['node_curve_xmaxymax']
                        for n in zip(nru,nlu):
                            cubit.cmd('equivalence node '+' '.join(str(x) for x in n)+' tolerance '+str(tol))
        #
        #
        cmd='topology check coincident node face all tolerance '+str(tol*2)+' nodraw brief result group "checkcoinc"' 
        cubit.silent_cmd(cmd)
        group_id_1=cubit.get_id_from_name("checkcoinc")
        if group_id_1 != 0:
            n1=cubit.get_group_nodes(group_id_1)
            if len(n1) != 0:
                print 'error, coincident nodes after the equivalence node command, check the tolerance'
                import sys
                sys.exit()
        cmd='group "negativejac" add quality hex all Jacobian high'
        cubit.cmd(cmd) 
        group_id_1=cubit.get_id_from_name("negativejac")
        n1=cubit.get_group_nodes(group_id_1)
        if len(n1) != 0:
            print 'error, negative jacobian after the equivalence node command, check the mesh'
    elif ckbound_method1 and  ckbound_method2 and len(filenames) != 1:
        block_list=cubit.get_block_id_list()
        i=-1
        for block in block_list:
            ty=cubit.get_block_element_type(block)
            if ty == 'HEX8':
                i=i+1
                if isinstance(merge_tolerance,list):
                    try:
                        tol=merge_tolerance[i]
                    except:
                        tol=merge_tolerance[-1]
                elif merge_tolerance:
                    tol=merge_tolerance
                else:
                    tol=1
                cmd='topology check coincident node face in hex in block '+str(block)+' tolerance '+str(tol)+' nodraw brief result group "b'+str(block)+'"'
                cubit.cmd(cmd)
                print cmd
                cmd='equivalence node in group b'+str(block)+' tolerance '+str(tol)
                cubit.cmd(cmd)
                print cmd
        if isinstance(merge_tolerance,list):
            tol=max(merge_tolerance)
        elif merge_tolerance:
            tol=merge_tolerance
        else:
            tol=1
        #
        #
        cmd='topology check coincident node face all tolerance '+str(tol)+' nodraw brief result group "checkcoinc"' 
        cubit.silent_cmd(cmd)
        group_id_1=cubit.get_id_from_name("checkcoinc")
        if group_id_1 != 0:
            n1=cubit.get_group_nodes(group_id_1)
            if len(n1) != 0:
                print 'error, coincident nodes after the equivalence node command, check the tolerance'
                import sys
                sys.exit()
        cmd='group "negativejac" add quality hex all Jacobian high'
        cubit.silent_cmd(cmd) 
        group_id_1=cubit.get_id_from_name("negativejac")
        n1=cubit.get_group_nodes(group_id_1)
        if len(n1) != 0:
            print 'error, negative jacobian after the equivalence node command, use --merge instead of --equivalence'




def collect(cpuxmin=0,cpuxmax=1,cpuymin=0,cpuymax=1,cpux=1,cpuy=1,cubfiles=False,ckbound_method1=False,ckbound_method2=False,merge_tolerance=None,curverefining=False,outfilename='totalmesh_merged',qlog=False,export2SPECFEM3D=False,listblock=None,listflag=None,outdir='.',add_sea=False,decimate=False,cpml=False,cpml_size=False,top_absorbing=False,hex27=False):
    #
    cubit.cmd('set journal error off')
    cubit.cmd('set verbose error off')
    collecting_merging(cpuxmin,cpuxmax,cpuymin,cpuymax,cpux,cpuy,cubfiles=cubfiles,ckbound_method1=ckbound_method1,ckbound_method2=ckbound_method2,merge_tolerance=merge_tolerance,decimate=decimate)
    cubit.cmd('set journal error on')
    cubit.cmd('set verbose error on')
    #
    if curverefining:
        block=1001 #topography
        refine_closecurve(block,curverefining,acis=True)
    #
    #
    if add_sea:
        block=1001
        add_sea_layer(block=block)

    outdir2='/'.join(x for x in outfilename.split('/')[:-1])
    if outdir2 == '': 
        outdir2=outdir+'/'
    else:
        outdir2=outdir+'/'+outdir2+'/'    
    
    import os
    try:
        os.makedirs(outdir2)
    except OSError:
        pass
    
    cubit.cmd('compress all')
    command="export mesh '"+outdir2+outfilename+".e' block all overwrite xml '"+outdir2+outfilename+".xml'"
    cubit.cmd(command)
    f=open(outdir2+'blocks.dat','w')
    blocks=cubit.get_block_id_list()
    #
    for block in blocks:
        name=cubit.get_exodus_entity_name('block',block)
        element_count = cubit.get_exodus_element_count(block, "block")
        nattrib=cubit.get_block_attribute_count(block)
        attr=[cubit.get_block_attribute_value(block,x) for x in range(0,nattrib)]
        ty=cubit.get_block_element_type(block)
        f.write(str(block)+' ; '+name+' ; nattr '+str(nattrib)+' ; '+' '.join(str(x) for x in attr)+' ; '+ty+' '+str(element_count)+'\n')
    f.close()
    #
    #
    cubit.cmd('set info echo journ off')
    cmd='del group all'
    cubit.silent_cmd(cmd)
    cubit.cmd('set info echo journ on')
    #
    command = "save as '"+outdir2+outfilename+".cub' overwrite"
    cubit.cmd(command)
    #
    print 'end meshing'
    #
    #
    if qlog:
        print '\n\nQUALITY CHECK.... ***************\n\n'
        import quality_log
        tq=open(outdir2+outfilename+'.quality','w')
        max_skewness,min_length=quality_log.quality_log(tq)
    #
    #
    #
    if export2SPECFEM3D:
        e2SEM(files=False,listblock=listblock,listflag=listflag,outdir=outdir,cpml=cpml,cpml_size=cpml_size,top_absorbing=top_absorbing,hex27=hex27)
                                            
def e2SEM(files=False,listblock=None,listflag=None,outdir='.',cpml=False,cpml_size=False,top_absorbing=False,hex27=False):
    import glob
    if files:
        filenames=glob.glob(files)
        for f in filenames:
            print f
            extension=f.split('.')[-1]
            if extension == 'cub':
                cubit.cmd('open "'+f+'"')
            elif extension== 'e':
                cubit.cmd('import mesh "'+f+'" no_geom')
            else:
                print extension
    if listblock and listflag:
        pass
    else:
        listblock=[]
        listflag=[]
        block_list=list(cubit.get_block_id_list())
        for block in block_list:
            ty=cubit.get_block_element_type(block)
            if 'HEX' in ty:
                listblock.append(block)
                #listflag.append(block)
        listflag=range(1,len(block_list)+1)  
    #       
    for ib,iflag in zip(listblock,listflag):
        cubit.cmd("block "+str(ib)+" attribute count 1")
        cubit.cmd("block "+str(ib)+" attribute index 1 "+ str(iflag)            )
    #
    import cubit2specfem3d
    cubit2specfem3d.export2SPECFEM3D(outdir,cpml=cpml,cpml_size=cpml_size,top_absorbing=top_absorbing,hex27=hex27)

def invert_dict(d):
     inv = {}
     for k,v in d.iteritems():
         keys = inv.setdefault(v, [])
         keys.append(k)
     return inv
     
     
     
def prepare_equivalence(nodes1,nodes2):
    cubit.cmd('set info off')
    cubit.cmd('set echo off')
    cubit.cmd('set journal off')
    length={}
    for ns in zip(nodes1,nodes2):
        cmd='group "tmpn" add edge in node '+' '.join(str(n) for n in ns )
        cubit.cmd(cmd)
        ge=cubit.get_id_from_name("tmpn")
        e1=cubit.get_group_edges(ge)
        lengthmin=1e9
        for e in e1:
            lengthmin=min(lengthmin,cubit.get_mesh_edge_length(e))
        length[ns]=lengthmin*.5
        cubit.cmd('delete group '+str(ge))
    minvalue=min(length.values())
    maxvalue=max(length.values())
    print 'min lentgh: ',minvalue,'max lentgh: ',maxvalue
    nbin= int((maxvalue/minvalue)/2.)+1
    factor=(maxvalue-minvalue)/nbin
    dic_new={}
    for k in length.keys():
        dic_new[k]=int((length[k]-minvalue)/factor)
    inv_length=invert_dict(dic_new)
    print inv_length.keys(),factor,minvalue
    ks=inv_length.keys()
    ks.sort()
    for k in range(0,len(inv_length.keys())-1):
        inv_length[ks[k]]=inv_length[ks[k]]+inv_length[ks[k+1]]
    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    cubit.cmd('set journal on')
    return factor,minvalue,inv_length


def merge_node_ck(n1,n2):
    factor,minvalue,inv_length=prepare_equivalence(n1,n2)
    
    cubit.cmd('set info off')
    cubit.cmd('set echo off')
    cubit.cmd('set journal off')
    cubit.cmd('set error off')
    
    for k in inv_length.keys()[:-1]:
        if len(inv_length[k]) > 0:
            cmd='equivalence node '+' '.join(' '.join(str(n) for n in x) for x in inv_length[k])+' tolerance '+str(k*factor+minvalue/2.)
            cubit.cmd(cmd)
            print 'equivalence '+str(len(inv_length[k]))+' couples of nodes -  tolerance '+str(k*factor+minvalue/2.)
             

    cubit.cmd('group "checkmerge" add node '+' '.join(str(n) for n in n1)+' '+' '.join(str(n) for n in n2))
    idg=cubit.get_id_from_name('checkmerge')
    remainnodes=cubit.get_group_nodes(idg)
    print 'from '+str(len(n1)+len(n2))+' nodes -> '+str(len(remainnodes)) +' nodes'
    if len(n1) != len(remainnodes):
        print 'equivalence '+str(len(remainnodes))+' couples of nodes -  tolerance '+str(minvalue/2.)
        cubit.cmd('set info on')
        cubit.cmd('set echo on')
        cubit.cmd('set journal on')
        cmd='equivalence node in group '+str(idg)+' tolerance '+str(minvalue/2.)
        cubit.cmd(cmd)
        cmd='block 3000 node in group '+str(idg)
        cubit.cmd(cmd)
        
    if len(n1) != len(remainnodes):
        cubit.cmd('export mesh "error_merging.e" dimension 3 block all overwrite')
        cubit.cmd('save as "error_merging.cub" dimension 3 block all overwrite')
        print 'error merging '
        if False:
            import sys
            sys.exit(2)
    
    cubit.cmd('delete group checkmerge')
    cubit.cmd('delete block 3000')
    
    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    cubit.cmd('set journal on')




def merge_node(n1,n2):
    factor,minvalue,inv_length=prepare_equivalence(n1,n2)

    cubit.cmd('set info off')
    cubit.cmd('set echo off')
    cubit.cmd('set journal off')


    for k in inv_length.keys()[:-1]:
        if len(inv_length[k]) > 0:
            cmd='equivalence node '+' '.join(' '.join(str(n) for n in x) for x in inv_length[k])+' tolerance '+str(k*factor+minvalue/2.)
            cubit.cmd(cmd)
            print 'equivalence '+str(len(inv_length[k]))+' couples of nodes -  tolerance '+str(k*factor+minvalue/2.)

    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    cubit.cmd('set journal on')


    
    

def prepare_equivalence_4(nodes1,nodes2,nodes3,nodes4):
    cubit.cmd('set info off')
    cubit.cmd('set echo off')
    cubit.cmd('set journal off')
    length={}
    nodes=[nodes1,nodes2,nodes3,nodes4]
    check=map(len,nodes)
    checked_nodes=[]
    for ind,iflag in enumerate(check):
        if iflag:
            checked_nodes=checked_nodes+nodes[ind]
    
    cmd='group "tmpn" add edge in node '+' '.join(str(n) for n in checked_nodes )
    cubit.cmd(cmd)
    ge=cubit.get_id_from_name("tmpn")
    e1=cubit.get_group_edges(ge)
    lengthmin=1e9
    for e in e1:
        lengthmin=min(lengthmin,cubit.get_mesh_edge_length(e))
        length[e]=lengthmin*.5
    cubit.cmd('delete group '+str(ge))
    try:
        minvalue=min(length.values())
        maxvalue=max(length.values())
    except:
        try:
            print nodes
            print 'edges ', e1
        except:
            pass
        minvalue=10.
        maxvalue=2000.
    print 'min lentgh: ',minvalue,'max lentgh: ',maxvalue
    nbin= int((maxvalue/minvalue)/2.)+1
    factor=(maxvalue-minvalue)/nbin
    dic_new={}
    for k in length.keys():
        dic_new[k]=int((length[k]-minvalue)/factor)
    inv_length=invert_dict(dic_new)
    print inv_length.keys(),factor,minvalue
    ks=inv_length.keys()
    ks.sort()
    for k in range(0,len(inv_length.keys())-1):
        inv_length[ks[k]]=inv_length[ks[k]]+inv_length[ks[k+1]]
    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    cubit.cmd('set journal on')
    return factor,minvalue,inv_length

def ording_z(nodes):
    def get_z(node):
        x,y,z = cubit.get_nodal_coordinates(node)
        return z
    d = [(get_z(node), node) for node in nodes]
    d.sort()
    return [x[1] for x in d] 

def merge_node_4(n1,n2,n3,n4,newmethod=True):
    if newmethod:
        print "merge node 4 side"
        n1o=ording_z(n1)
        n2o=ording_z(n2)
        n3o=ording_z(n3)
        n4o=ording_z(n4)
        for ln in zip(n1o,n2o,n3o,n4o):
            cmd='equivalence node '+' '.join(str(n) for n in ln) +' tolerance 10000 '
            cubit.cmd(cmd)
        #    
        #allnodes=n1+n2+n3+n4
        #print allnodes
        #for n in allnodes:
        #    print n
        #cmd='equivalence node '+' '.join(str(n) for n in allnodes) +' tolerance 10 ' 
        #cubit.cmd(cmd)
    else:
        factor,minvalue,inv_length=prepare_equivalence_4(n1,n2,n3,n4)
        
        for k in inv_length.keys()[:-1]:
            if len(inv_length[k]) > 1:
                try:
                    for x in inv_length[k]:
                        if type(x) is not list:
                            x=[x]
                        else:
                            pass
                    cmd='equivalence node '+' '.join(' '.join(str(n) for n in x) )+' tolerance '+str(k*factor+minvalue/2.)
                except:
                    print k,"***************************************** s"
                    print inv_length[k]
                    
                cubit.cmd(cmd)
                print 'equivalence '+str(len(inv_length[k]))+' couples of nodes -  tolerance '+str(k*factor+minvalue/2.)
            if len(inv_length[k]) == 1:
                cmd='equivalence node '+' '.join(' '.join(str(n) for n in  inv_length[k]))+' tolerance '+str(k*factor+minvalue/2.)
                cubit.cmd(cmd)
                print 'equivalence '+str(len(inv_length[k]))+' couples of nodes -  tolerance '+str(k*factor+minvalue/2.)

