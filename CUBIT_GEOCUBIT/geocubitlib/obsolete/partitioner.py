#############################################################################
# partitioner.py                                                    
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

def partitioner_obsolete():
    """create the partitioner"""
    import initializing as initializing
    cubit                       = initializing.initializing_cubit()
    cfg                         = initializing.initializing_cfg()
    
    import sys, os
    from sets import Set
    from utilities import load_curves,project_curve
    
    ############
    #
    command = "reset"
    cubit.cmd(command)
    #
    if cfg.debug:
        command = "set info on"
        cubit.cmd(command)
        command = "set logging on file '"+cfg.working_dir+"/debug.log'"
        cubit.cmd(command)
    #
    xmin_surf=cfg.xmin
    xmax_surf=cfg.xmax
    ymin_surf=cfg.ymin
    ymax_surf=cfg.ymax
    #
    vstart=cubit.get_last_id("vertex")+1
    command = "create vertex "+str(xmin_surf)+" "+str(ymin_surf)+" "+str(cfg.top_partitioner)
    cubit.cmd(command)
    command = "create vertex "+str(xmin_surf)+" "+str(ymax_surf)+" "+str(cfg.top_partitioner)
    cubit.cmd(command)
    command = "create vertex "+str(xmax_surf)+" "+str(ymax_surf)+" "+str(cfg.top_partitioner)
    cubit.cmd(command)
    command = "create vertex "+str(xmax_surf)+" "+str(ymin_surf)+" "+str(cfg.top_partitioner)
    cubit.cmd(command)
    vend=cubit.get_last_id("vertex")
    #
    command = "create surface vertex " +str(vstart)+" to "+str(vend)
    cubit.cmd(command)
    #
    top_surface=cubit.get_last_id("surface")
    top_volume=cubit.get_last_id("volume")
    #
    external_curve=cubit.get_relatives("surface", top_surface, "curve")
    external_curve=str(list(external_curve)).replace("["," ").replace("]"," ")
    print external_curve
    print cfg.outline_curve,cfg.transition_curve,cfg.partitioning_curve,cfg.internal_curve
    #
    try:
        outline_curve=load_curves(cfg.outline_curve,top_surface)
    except:
        outline_curve=None
    try:
        transition_curve=load_curves(cfg.transition_curve,top_surface)
    except:
        transition_curve=None
    try:
        internal_curve=load_curves(cfg.internal_curve,top_surface)
    except:
        internal_curve=None
    try:
        partitioning_curve=load_curves(cfg.partitioning_curve,top_surface)
    except:
        partitioning_curve=None
    print outline_curve,transition_curve,partitioning_curve,internal_curve
    #
    command = "save as '"+ cfg.working_dir+ "/partitioner_debug.cub' overwrite"
    cubit.cmd(command)    
    #
    #
    if internal_curve:
        print 'internal'
        tmp_curve=cubit.get_last_id("curve")
        cubit.cmd('imprint volume all with curve '+str(internal_curve))
        tmp_curve_after=cubit.get_last_id("curve")
        internal_curve=str(range(tmp_curve+1,tmp_curve_after+1)).replace("["," ").replace("]"," ")
    #
    if transition_curve:
        tmp_curve=cubit.get_last_id("curve")
        cubit.cmd('imprint volume all with curve '+str(transition_curve))
        tmp_curve_after=cubit.get_last_id("curve")
        transition_curve=str(range(tmp_curve+1,tmp_curve_after+1)).replace("["," ").replace("]"," ")
    #
    if partitioning_curve:
        cubit.cmd('imprint volume all with curve '+str(partitioning_curve)) 
    #
    cubit.cmd('delete curve all')
    cubit.cmd('imprint all')
    cubit.cmd('merge surface all')
    cubit.cmd('merge curve all')
    #
    mer=cubit.get_error_count()
    command = "save as '"+ cfg.working_dir+ "/partitioner.cub' overwrite"
    cubit.cmd(command)    
    #
    ls=cubit.parse_cubit_list("surface","all")
    print ls
    s_out=Set(ls)
    #
    s_t=Set()
    if transition_curve:
        for s in transition_curve.split(','):
            s=int(s)
            #print Set((cubit.get_relatives("curve",s,"surface")))
            s_t=s_t | Set((cubit.get_relatives("curve",s,"surface")))
    #
    print 't',s_t
    s_inside=Set()
    if internal_curve:
        for s in internal_curve.split(','):
            s=int(s)
            print s,Set((cubit.get_relatives("curve",s,"surface")))
            s_inside=s_inside | Set((cubit.get_relatives("curve",s,"surface")))
    #
    print 'inside', s_inside
    if not internal_curve and transition_curve:
        s_out = s_out - s_t
        s_t = s_t - s_out
    else:
        s_in  = s_inside - s_t
        s_t   = s_inside - s_in
        s_out = s_out - s_t
        s_out = s_out - s_in
    #
    surf_out=list(s_out)
    surf_out_txt=str(surf_out).replace("["," ").replace("]"," ")
    surf_transition=list(s_t)
    surf_transition_txt=str(surf_transition).replace("["," ").replace("]"," ")
    surf_in=list(s_in)
    surf_in_txt=str(surf_in).replace("["," ").replace("]"," ")
    print 'in ',surf_in
    print 'out ',surf_out
    print 'transition ',surf_transition
    #
    cubit.cmd('surface '+surf_in_txt+' size '+str(cfg.size_small))
    cubit.cmd('surface '+str(surf_in_txt)+' scheme pave')
    cubit.cmd('mesh surface '+str(surf_in_txt))
    #
    #cubit.cmd('curve '+transition_curve+' size '+str(cfg.size_large))
    cubit.cmd('curve '+transition_curve+' interval 1')
    command = "mesh curve "+transition_curve
    cubit.cmd(command)
    cubit.cmd('surface '+str(surf_out_txt)+' size '+str(cfg.size_large))
    cubit.cmd('surface '+str(surf_out_txt)+' scheme pave')
    cubit.cmd('mesh surface '+str(surf_out_txt))
    #
    #TODO: make automatic
    if  cfg.manual_adj:
        command = "save as '"+cfg.working_dir+'/'+cfg.file_manual_adjustment+"' overwrite"
        cubit.cmd(command)
        print '*************************'
        print command
        print 'manual node adjustment'
        q=0
        for l in ls:
            q=q+cubit.get_surface_element_count(l)
        print q, ' quads ---> CPUs'
        sys.exit()
    elif cfg.play_adj:
        print cfg.play_adj
        for command in cfg.play_adj.split('\n'):
            cubit.cmd(command)
            print command
    elif cfg.no_adj:
        print 'no adjustements'
    else:
        print 'please check the configuration: one of manual_adj,play_adj,no_adjust should be True'
    
    
    #
    command = "group 'list_edge_out' add edge in face in surf "+str(surf_out_txt)
    cubit.cmd(command)
    command = "group 'list_edge_t' add edge in face in surf "+str(surf_transition_txt)
    cubit.cmd(command)
    try:
        group=cubit.get_id_from_name("list_edge_in")
    except:
        command = "group 'list_edge_t' add edge in surf "+str(surf_transition_txt)
        cubit.cmd(command)
        group=cubit.get_id_from_name("list_edge_in")
    list_edge_in=cubit.get_group_edges(group)
    group=cubit.get_id_from_name("list_edge_out")
    list_edge_out=cubit.get_group_edges(group)
    group=cubit.get_id_from_name("list_edge_t")
    list_edge_t=cubit.get_group_edges(group)   
    edge_in=Set(list_edge_in)
    edge_out=Set(list_edge_out)
    edge_transition=Set(list_edge_t)
    edge_in=edge_in-edge_transition
    edge_out=edge_out-edge_transition
    #
    last_curve_start=cubit.get_last_id("curve")+1
    for edge in edge_in:
        nodes=cubit.get_connectivity("Edge",edge)
        command = "create curve spline node "+str(list(nodes))
        command = command.replace("["," ").replace("]"," ")
        cubit.cmd(command)
        last_curve=cubit.get_last_id("curve")
        command = "group 'l_curve_in' add curve "+str(last_curve)
        cubit.cmd(command)
    #
    for edge in list_edge_out:
        nodes=cubit.get_connectivity("Edge",edge)
        command = "create curve spline node "+str(list(nodes))
        command = command.replace("["," ").replace("]"," ")
        cubit.cmd(command)
        last_curve=cubit.get_last_id("curve")
        command = "group 'l_curve_out' add curve "+str(last_curve)
        cubit.cmd(command)
    #
    cubit.cmd('reset vol all')
    #
    igroup=cubit.get_id_from_name("l_curve_in")
    l_curve_in=cubit.get_group_curves(igroup)
    command = "delete group "+ str(igroup)
    cubit.cmd(command)
    igroup=cubit.get_id_from_name("l_curve_out")
    l_curve_out=cubit.get_group_curves(igroup)
    command = "delete group "+ str(igroup)
    cubit.cmd(command)
    #
    last_surface_before=cubit.get_last_id("surface")
    command = "imprint vol "+str(top_volume)+" with curve "+str(list(l_curve_in))
    command = command.replace("["," ").replace("]"," ")
    cubit.cmd(command)
    last_surface_after=cubit.get_last_id("surface")
    command = "group 'list_surface_in' add surface "+surf_in_txt+" "+str(last_surface_before+1)+" to "+str(last_surface_after)
    cubit.cmd(command)
    #
    last_surface_before=cubit.get_last_id("surface")                                                                            
    command = "imprint vol "+str(top_volume)+" with curve "+str(list(l_curve_out))                                              
    command = command.replace("["," ").replace("]"," ")
    cubit.cmd(command)
    last_surface_after=cubit.get_last_id("surface")
    command = "group 'list_surface_out' add surface "+surf_out_txt+" "+str(last_surface_before+1)+" to "+str(last_surface_after)
    cubit.cmd(command)
    #
    command = "group 'list_surface_transition' add surface "+surf_transition_txt
    cubit.cmd(command)
    #
    command = "del curve all"
    cubit.cmd(command)
    #
    cubit.cmd('merge surf in volume '+str(top_volume))
    cubit.cmd('list group all')
    command = "compress all"
    cubit.cmd(command)
    #TODO filename general
    command = "save as '"+cfg.working_dir+'/'+cfg.file_partitioner_cub+"' overwrite"
    cubit.cmd(command)
    print command
    cubit.print_surface_summary_stats()
    print 'saved'
    print
    ls=cubit.parse_cubit_list("surface","all")
    print '@@@@@@@@@@@@@@@@@@@@@@ ', ls, ' CORES requested ' 
    

def partitioner_map_obsolete(partitioner_filename=None):
    """return the map of the partioner surfaces associated to the cpu, and exclude some surface

    usage: list_surf_all,list_surface_in,list_surface_out,list_surface_transition,list_surface_nomesh,pmap=partitioner_map()

    the map is a dictionary of surfmap class objects:

        obj=surfmap(iproc,list_surface_all)

        obj.iproc = mpi processor
        obj.idsurf= surface id (the surface is saved in a ACIS (.sat) file called "surf_{idsurf}.sat")                                      
        obj.s_adj = get the id list of the surfaces adjacent to obj.idsurf
        obj.p_adj = get the id list of the processor associated to the surface adjacent to obj.idsurf
        obj.point = get the coordinate of the center_point of the curves shared by the adjacent surfaces
        obj.all   = list of all the previous obj properties [self.iproc,self.idsurf,self.s_adj,self.p_adj,self.point]


    for example:

        iproc processor -> id surf from list_surf_all / [list id surf adjacent] [list id proc adjacent]: 
        0   -> 4/[5, 1, 15] [1, 22, 11]
        1   -> 5/[4, 16, 17, 1] [0, 12, 13, 22]
        2   -> 6/[19, 22, 21] [15, 18, 17]
        3   -> 7/[8, 14, 9] [4, 10, 5]
        .....

    """
    import initializing as initializing
    cubit                       = initializing.initializing_cubit()
    cfg                         = initializing.initializing_cfg()
    #
    mpiflag,iproc,numproc,mpi   = initializing.initializing_mpi()
    from mpi_geocubit import mpiprint
    #
    #
    import sys, os
    from sets import Set
    #
    def get_surf_adj(iproc,list_surf_all):
          import cubit
          l=[]
          l=list(cubit.get_adjacent_surfaces("surface", list_surf_all[iproc]))
          l.remove( list_surf_all[iproc])
          return l
    #
    def get_proc_adj(iproc,list_surf_all,list_surf_adj):
          l=[]
          for s in list_surf_adj:
                l.append(list_surf_all.index(s))
          return l
    #
    def get_x(idsurf,list_surf_adj):
          import cubit
          from sets import Set
          l=[]
          for s in list_surf_adj:
                c=Set(cubit.get_relatives("surface",s,"curve")) & Set(cubit.get_relatives("surface",idsurf,"curve"))
                p=cubit.get_center_point("curve",list(c)[0])
                l.append(p)
          return l
    #
    #
    class surfmap:
          def __init__(self,iproc,list_surf_all):
              self.iproc=iproc
              self.idsurf=list_surf_all[iproc]
              self.s_adj=get_surf_adj(iproc,list_surf_all)
              self.p_adj=get_proc_adj(iproc,list_surf_all,self.s_adj)
              self.point=get_x(self.idsurf,self.s_adj) 
              self.all=[self.iproc,self.idsurf,self.s_adj,self.p_adj,self.point]     
          def __repr__(self):
              msg=str(list(self.all))
              return msg
    #
    command = "reset"
    cubit.cmd(command)
    #
    if cfg.create_partitioner:
       partitioner()
    elif os.path.exists(cfg.file_partitioner_cub):
       command = "import cubit '"+cfg.file_partitioner_cub+"'"
       cubit.cmd(command)
    elif os.path.exists(cfg.working_dir+'/'+cfg.file_partitioner_cub):
       command = "import cubit '"+cfg.working_dir+"/"+cfg.file_partitioner_cub+"'"
       cubit.cmd(command)
       
       
       
       
       
    #
    #
    #command = "imprint curve all"
    #cubit.cmd(command)
    #command = "imprint surface all"
    #cubit.cmd(command)
    #command = "merge surf all"
    #cubit.cmd(command)
    #command = "imprint all"
    #cubit.cmd(command)
    #command = "merge all"
    #cubit.cmd(command)
    command = "compress all"
    cubit.cmd(command)
    #
    group=cubit.get_id_from_name("list_surface_in")
    list_surface_in=cubit.get_group_surfaces(group)
    group=cubit.get_id_from_name("list_surface_out")
    list_surface_out=cubit.get_group_surfaces(group)
    group=cubit.get_id_from_name("list_surface_transition")
    list_surface_transition=cubit.get_group_surfaces(group)
    group=cubit.get_id_from_name("list_surface_nomesh")
    list_surface_nomesh=cubit.get_group_surfaces(group)
    cubit.cmd('list group all')
    print list_surface_in,list_surface_transition,list_surface_out,list_surface_nomesh
    
    #
    list_surface_all=list(list_surface_out+list_surface_in+list_surface_transition+list_surface_nomesh)
    list_surface_decomposer=list(list_surface_out+list_surface_in)
    #
    if os.path.exists(cfg.working_dir) and (iproc == 0 or cfg.single):
        command = "export acis '"+cfg.working_dir+"/surf_all.sat' surface all overwrite"
        cubit.cmd(command)
        #cubit_error_stop(iproc,command,ner)
        for i in list_surface_all:
            command = "export acis '"+cfg.working_dir+"/surf_"+str(i)+".sat' surface "+str(i)+" overwrite"
            cubit.cmd(command)
            #cubit_error_stop(iproc,command,ner)
    #
    pmap={}
    for i in range(0,len(list_surface_all)):
        pmap[i]= surfmap(i,list_surface_all)
    #
    #
    text='***********************\nthe partioner is composed of ***'+str(len(list_surface_all))+ '*** surfaces\nit requires ***'+str(len(list_surface_all))+'*** cpus\n******************'
    mpiprint(mpiflag,text)
    #
    if not mpiflag:
        print 'ATTENTION***********************************************'
        print 'id proc = ',iproc,' -> surf ',list_surface_all[iproc]
        print 'only the volume correspondent to surf_',str(list_surface_all[iproc]),'.sat will be mesh'
        print 'processor -> list_surf_all / parallel map: '
        for i in range(0,len(list_surface_all)):
                flag=''
                if list_surface_all[i] in list_surface_in:
                    flag='( in  )'
                elif list_surface_all[i] in list_surface_out:
                    flag='( out )'
                elif list_surface_all[i] in list_surface_transition:
                    flag='( t   )'
                elif list_surface_all[i] in list_surface_nomesh:
                    flag='( no   )'
                print i,' ',flag,' -> ', list_surface_all[i],'/', surfmap(i,list_surface_all).s_adj, surfmap(i,list_surface_all).p_adj
        print '********************************************************'
    else:
        if iproc == 0:
            partxt=open(cfg.working_dir+"/partitioner_map.txt","w")
            partxt.write('processor -> list_surf_all / parallel map: \n')
            for i in range(0,len(list_surface_all)):
                flag=''
                if list_surface_all[i] in list_surface_in:
                    flag='( in  )'
                if list_surface_all[i] in list_surface_out:
                    flag='( out )'
                if list_surface_all[i] in list_surface_transition:
                    flag='( t   )'
                if list_surface_all[i] in list_surface_nomesh:
                    flag='( no   )'
                partxt.write(str(i)+' '+flag+' -> '+ str(list_surface_all[i])+'/'+ str(surfmap(i,list_surface_all).s_adj)+' '+str(surfmap(i,list_surface_all).p_adj)+'\n')
            partxt.close()
    #
    #
    text='***********************\nthe partioner is composed of ***'+str(len(list_surface_all))+ '*** surfaces\nit requires ***'+str(len(list_surface_all))+'*** cpus\n******************'
    mpiprint(mpiflag,text)
    #
    #TODO make recursive so I can mesh with a lower number of cpu
    if numproc < len(list_surface_all) and mpiflag:
        text='ATTENTION**********************************************'
        mpiprint(mpiflag,text)
        for i in range(numproc+1,len(list_surface_all)):
            text='the volume correspondent to surf_'+str(list_surface_all[iproc])+'.sat will ***NOT*** be mesh'
            mpiprint(mpiflag,text)
        text='********************************************************'
        mpiprint(mpiflag,text)
    if numproc > len(list_surface_all):
        text='number of processors greater than the number of slices....'
        mpiprint(mpiflag,text)
        raise NameError, 'error'    
    #
    #
    return list_surface_all,list_surface_in,list_surface_out,list_surface_transition,list_surface_nomesh,pmap




def partitioning_volume_obsolete():
    import sys,os
    import menu as menu
    import initializing as initializing
    import partitioner as partitioner
    from sets import Set
    #
    mpiflag,iproc,numproc,mpi   = initializing.initializing_mpi()
    from mpi_geocubit import mpiprint
    #
    numpy                       = initializing.initializing_numpy()
    cubit                       = initializing.initializing_cubit()
    cfg                         = initializing.initializing_cfg()

    def sortf(x,y):
        return cmp(x[1],y[1])
    #
    #
    list_surface_all,list_surface_in,list_surface_out,list_surface_transition,list_surface_nomesh,pmap=partitioner.partitioner_map(cfg.file_partitioner_cub)
    #
    #
    command = "reset"
    cubit.cmd(command)
    #
    id_partitioner_surf=list_surface_all[iproc]
    command = "import acis '"+cfg.working_dir+"/surf_"+str(id_partitioner_surf)+".sat'"
    cubit.cmd(command)
    last_surface=cubit.get_last_id("surface")
    #
    top_box=cubit.get_bounding_box("surface",last_surface)
    s_xmin=cubit.get_bounding_box("surface",last_surface)[0]
    s_xmax=cubit.get_bounding_box("surface",last_surface)[1]
    s_ymin=cubit.get_bounding_box("surface",last_surface)[3]
    s_ymax=cubit.get_bounding_box("surface",last_surface)[4]
    command = "del surf "+str(last_surface)
    cubit.cmd(command)
    #
    if cfg.volume_file:
        command = "import cubit '"+cfg.volume_file+"'"
        cubit.cmd(command)
    else:
        import local_volume
        if iproc == 0:
            coordx_0,coordy_0,elev_0,nx_0,ny_0=local_volume.read_grid()
            coordx=mpi.bcast(coordx_0)
            coordy=mpi.bcast(coordy_0)
            elev=mpi.bcast(elev_0)
            nx=mpi.bcast(nx_0)
            ny=mpi.bcast(ny_0)
        else:
            coordx=mpi.bcast()
            coordy=mpi.bcast()
            elev=mpi.bcast()
            nx=mpi.bcast()
            ny=mpi.bcast()
        extract_volume(s_xmin,s_ymin,s_xmax,s_ymax,coordx,coordy,elev,nx,ny)
    #cubit_error_stop(iproc,command,ner)
    #
    list_vol=cubit.parse_cubit_list("volume","all")
    command = "import acis '"+cfg.working_dir+"/surf_"+str(id_partitioner_surf)+".sat'"
    cubit.cmd(command)
    last_surface=cubit.get_last_id("surface")
    #
    init_n_vol=len(list_vol)
    zmin_box=cubit.get_total_bounding_box("volume",list_vol)[6] #it is the z_min of the box ... box= xmin,xmax,d,ymin,ymax,d,zmin...    
    xmin_box=cubit.get_total_bounding_box("volume",list_vol)[0]
    xmax_box=cubit.get_total_bounding_box("volume",list_vol)[1]
    ymin_box=cubit.get_total_bounding_box("volume",list_vol)[3]
    ymax_box=cubit.get_total_bounding_box("volume",list_vol)[4]
    zmax_box=cubit.get_bounding_box("surface",last_surface)[7]
    #
    distance=2.*(zmax_box-zmin_box)
    #
    command = "composite create curve in surf "+str(last_surface)
    cubit.cmd(command)
    ner=cubit.get_error_count()
    command = "webcut volume "+str(list(list_vol))+" sweep surface "+str(last_surface)+" vector 0 0 -1 distance "+str(distance)
    command = command.replace("["," ").replace("]"," ")
    cubit.cmd(command)
    #cubit_error_stop(iproc,command,ner)
    #
    #
    #delete the vol 
    list_vol=cubit.parse_cubit_list("volume","all")
    set_vol=Set([])
    dbox=[]
    for id_vol in list_vol:
        lsurf=cubit.get_relatives("volume",id_vol,"surface")
        set_vol.add(len(lsurf))
        d=cubit.get_bounding_box("volume",id_vol)[1]-cubit.get_bounding_box("volume",id_vol)[0]
        dbox.append([id_vol,d])
    dbox.sort(sortf)
    #TODO explain....
    if id_partitioner_surf in list_surface_transition:
        if len(set_vol) != 3: #in the transition there is a hole so there are 3 classes of volume 
           txt='transition - imprint of partitioner surface not consistent'+str(id_partitioner_surf)
           command = "comment '"+txt+"'"
           cubit.cmd(command)
           raise NameError, txt
        for obj in dbox[0:init_n_vol]:
            command = "delete vol "+str(obj[0])
            cubit.cmd(command)
        for obj in dbox[3*init_n_vol-init_n_vol:3*init_n_vol]:
            command = "delete vol "+str(obj[0])
            cubit.cmd(command)
    else:
        if len(set_vol) != 2:
           txt='in-out-nomesh - imprint of partitioner surface not consistent'+str(id_partitioner_surf)
           command = "comment '"+txt+"'"
           cubit.cmd(command)
           raise NameError, txt 
        min_surf=min(set_vol)
        for id_vol in list_vol:
            lsurf=cubit.get_relatives("volume",id_vol,"surface")
            if  len(lsurf) != min_surf: 
                command = "delete vol "+str(id_vol) # I'm deleting the volume which doesn't have the minimum number of surfaces
                cubit.cmd(command)
                #cubit_error_stop(iproc,command,ner)
    #
    command = "delete surf "+str(last_surface)
    cubit.cmd(command)

    tol=cfg.size_small/10.
    command = "merge tolerance "+str(tol)
    cubit.cmd(command); print command
    command = "imprint tolerant volume all"
    cubit.cmd(command)
    #cubit_error_stop(iproc,command,ner)
    #command = "imprint all"
    #cubit.cmd(command)
    ##cubit_error_stop(iproc,command,ner)
    #command = "merge all"
    #cubit.cmd(command)
    #cubit_error_stop(iproc,command,ner)    
    command = "merge surf all"
    cubit.cmd(command)
    #cubit_error_stop(iproc,command,ner)
    command = "compress vol all"
    cubit.cmd(command)
    #cubit_error_stop(iproc,command,ner)
    command = "compress surf all"
    cubit.cmd(command)
    #cubit_error_stop(iproc,command,ner)
    if id_partitioner_surf in list_surface_out:
       command = "save as '"+cfg.working_dir+"/"+"out_vol_"+str(iproc)+".cub' overwrite"
       cubit.cmd(command)    
    elif id_partitioner_surf in list_surface_in:
        command = "save as '"+cfg.working_dir+"/"+"in_vol_"+str(iproc)+".cub' overwrite"
        cubit.cmd(command)
        #cubit_error_stop(iproc,command,ner)
    elif id_partitioner_surf in list_surface_transition:
        command = "save as '"+cfg.working_dir+"/"+"transition_vol_"+str(iproc)+".cub' overwrite"
        cubit.cmd(command)             
    elif id_partitioner_surf in list_surface_transition:
        command = "save as '"+cfg.working_dir+"/"+"nomesh_vol_"+str(iproc)+".cub' overwrite"
        cubit.cmd(command)                            
    else:
        raise NameError, 'surface not in list'    


def refinement_basin_obsolete(iproc,id_partitioner_surf,list_surface_out,list_surface_transition,top_surface,top_surface_add):
    import initializing as initializing
    cubit                       = initializing.initializing_cubit()
    cfg                         = initializing.initializing_cfg()
    if  id_partitioner_surf not in list_surface_out:
        from select_in_out import select_in_out
        cubit.cmd('save as "'+cfg.working_dir+'/view_before_basinrefinement_'+str(iproc)+'.cub" overwrite ')
        command = "set node constraint off"
        cubit.cmd(command)
        vc=cubit.get_volume_element_count(vol[-1].ID)
        sc=cubit.get_surface_element_count(top_surface)
        if cfg.cut_outline: 
            sc2=cubit.get_surface_element_count(top_surface_add)
        else:
            sc2=0
        #
        ntime=int(vc/(sc+sc2))-2
        #
        if not cfg.cut_outline and id_partitioner_surf in list_surface_transition:
            command = "import acis '"+cfg.outline_curve+"'"
            cubit.cmd(command)
            #cubit_error_stop(iproc,command,ner)
            curve_outline=cubit.get_last_id("curve")
            node_inside=select_in_out(top_surface,curve_outline)#!!!!!
            command = "group 'face_selected' add face in node "+str(list(node_inside)).replace("["," ").replace("]"," ")
            cubit.cmd(command)
            id_group=cubit.get_id_from_name('face_selected')
            command = "group "+str(id_group)+"remove face in (surf all except surf "+str(top_surface)+")"
            cubit.cmd(command)
        else:
            command = "group 'face_selected' add face in surf "+str(top_surface)
            cubit.cmd(command)
        #
        id_group=cubit.get_id_from_name('face_selected')
        if ntime > 0:
            command = "group 'hex_refining' add hex propagate face in group "+str(list(id_group))+"times "+str(ntime)
            command = command.replace("["," ").replace("]"," ")
            cubit.cmd(command)
            id_group=cubit.get_id_from_name('hex_refining')
            cubit.cmd('block 1 face in surf all ')
            cubit.cmd('set large exodus file off')
            cubit.cmd('save as "'+cfg.working_dir+'/before_basinref_view_side_'+str(iproc)+'.cub" overwrite ')
            command = "refine hex in group "+str(id_group)+'  numsplit 1 bias 1.0 depth 1 no_smooth'
            cubit.cmd(command)
        else:
            command = "refine face in group "+str(id_group)+'  numsplit 1 bias 1.0 depth 1 no_smooth'
            cubit.cmd(command)



def interpolate_basin_obsolete(top_surface,top_surface_add):
    import initializing as initializing
    cubit                       = initializing.initializing_cubit()
    cfg                         = initializing.initializing_cfg()
    #
    from select_in_out import select_in_out,select_in_out_global
    #
    id_group=cubit.get_id_from_name('hex_refining')
    command = "del group "+str(id_group)
    cubit.cmd(command)
    id_group=cubit.get_id_from_name('face_selected')
    command = "del group "+str(id_group)
    cubit.cmd(command)
    #
    vc=cubit.get_volume_element_count(vol[-1].ID)
    sc=cubit.get_surface_element_count(top_surface)
    if cfg.cut_outline: 
        sc2=cubit.get_surface_element_count(top_surface_add)
    else:
        sc2=0
    #
    ntime=int(vc/(sc+sc2))-2
    command = "import acis '"+cfg.outline_curve+"'"
    cubit.cmd(command)
    curve_outline=cubit.get_last_id("curve")
    node_inside=select_in_out(top_surface,curve_outline)
    command = "group 'face_selected' add face in node "+str(list(node_inside)).replace("["," ").replace("]"," ")
    cubit.cmd(command)
    id_group=cubit.get_id_from_name('face_selected')
    command = "group "+str(id_group)+"remove face in (surf all except surf "+str(top_surface)+")"
    cubit.cmd(command)
    if ntime > 0:
        command = "group 'hex_refining' add hex propagate face in group "+str(list(id_group))+"times "+str(ntime)
        command = command.replace("["," ").replace("]"," ")
        cubit.cmd(command)
    else:
        command = "group 'hex_refining' add hex in face in group "+str(list(id_group))
        command = command.replace("["," ").replace("]"," ")
        cubit.cmd(command)
    #
    id_group=cubit.get_id_from_name('hex_refining')
    list_hex=cubit.get_group_hexes(id_group)
    #
    command = "group 'list_node_inbasinhex' add node in hex "+str(list(list_hex))
    command = command.replace("["," ").replace("]"," ") 
    cubit.cmd(command)
    group=cubit.get_id_from_name("list_node_inbasinhex")
    list_node_inbasinhex=cubit.get_group_nodes(group)
    command = "delete group "+ str(group)
    cubit.cmd(command)
    #
    command = "import acis '"+cfg.basin_surf_acis_file+"'"
    cubit.cmd(command)
    basin_surface=cubit.get_last_id("surface")
    #
    nodes_inside_basin=select_in_out_global(list_node_inbasinhex,basin_surface)
    command = "nodeset 1 "+str(nodes_inside_basin)
    command = command.replace("["," ").replace("]"," ")
    cubit.cmd(command)







def mesh_partitioner_obsolete():
        import sys,os
        import menu as menu
        import initializing as initializing
        import partitioner as partitioner
        #
        mpiflag,iproc,numproc,mpi   = initializing.initializing_mpi()
        from mpi_geocubit import mpiprint
        #
        cubit                       = initializing.initializing_cubit()
        cfg                         = initializing.initializing_cfg()
        class cubitvolume:
              def __init__(self,ID,intervalv,centerpoint,dimension):
                  self.ID=ID
                  self.intervalv=intervalv
                  self.centerpoint=centerpoint
                  self.dim=dimension

              def __repr__(self):
                  msg="(vol:%3i, vertical interval: %4i, centerpoint: %8.2f)" % (self.ID, self.intervalv,self.centerpoint)
                  return msg       

        def by_z(x,y):
            return cmp(x.centerpoint,y.centerpoint)
        #
        list_surface_all,list_surface_in,list_surface_out,list_surface_transition,list_surface_nomesh,pmap=partitioner.partitioner_map(cfg.file_partitioner_cub)

        #
        command = "reset"
        cubit.cmd(command)
        #
        id_partitioner_surf=list_surface_all[iproc]
        #

        if id_partitioner_surf in list_surface_out:
           command = "import cubit '"+cfg.working_dir+"/out_vol_"+str(iproc)+".cub' "
           cubit.cmd(command)    
        elif id_partitioner_surf in list_surface_in:
            command = "import cubit '"+cfg.working_dir+"/in_vol_"+str(iproc)+".cub' "
            cubit.cmd(command)
        elif id_partitioner_surf in list_surface_transition:
            command = "import cubit '"+cfg.working_dir+"/transition_vol_"+str(iproc)+".cub' "
            cubit.cmd(command)
        elif id_partitioner_surf in list_surface_nomesh:
            command = "import cubit '"+cfg.working_dir+"/nomesh_vol_"+str(iproc)+".cub' "
            cubit.cmd(command)                                     
        else:
            raise NameError, 'surface not in list'
        if iproc == 0: print 'import...',command
        #
        list_vol=cubit.parse_cubit_list("volume","all")
        nvol=len(list_vol)                                 
        vol=[]
        for id_vol in list_vol:
            p=cubit.get_center_point("volume",id_vol)
            vol.append(cubitvolume(id_vol,1,p[2],0))
        vol.sort(by_z)
        for id_vol in range(0,nvol):
            vol[id_vol].intervalv=cfg.iv_interval[id_vol]
        #
        #TRANSITION SURF -     
        if  cfg.cut_outline and id_partitioner_surf in list_surface_transition:
            text='outline cutting.......'
            #mpiprint(mpiflag,text)
            command = "comment '"+"outline cutting......."+"'"
            cubit.cmd(command)
            command = "import acis '"+cfg.outline_curve+"'"
            cubit.cmd(command)
            #cubit_error_stop(iproc,command,ner)
            outline_curve=cubit.get_last_id("curve")
            #
            #
            lsurf=cubit.get_relatives("volume",vol[-1].ID,"surface")
            normal=[]
            tres=0.2
            for k in lsurf:
                normal=cubit.get_surface_normal(k)
                if abs(normal[2] - 1) <= tres:
                    top_surf=k
            #
            command = "project curve "+str(outline_curve)+ "onto surface " +str(top_surf)
            cubit.cmd(command)
            #
            outline_curve=cubit.get_last_id("curve")
            #
            command = "imprint vol "+str(vol[nvol-1].ID)+" with curve "+str(outline_curve)
            cubit.cmd(command)
            #
            command = "delete curve all"
            cubit.cmd(command)
            command = "imprint all"
            cubit.cmd(command)
            #cubit_error_stop(iproc,command,ner)
            command = "merge all"
            cubit.cmd(command)
        #
        #
        surf_vertical=[]
        surf_or=[]
        top_surface=0
        top_surface_add=''
        bottom_surface=0
        for l in range(0,nvol):
            id_vol=vol[l].ID
            lsurf=cubit.get_relatives("volume",id_vol,"surface")
            znormal=[]
            tres=cfg.tres
            icheck=0
            list_curve_vertical=[]
            list_curve_or=[]
            surf_vertical_tmp=[]
            source_sweep=0
            target_sweep=0
            for k in lsurf:
                normal=cubit.get_surface_normal(k)
                center_point = cubit.get_center_point("surface", k)
                if normal[2] >= -1*tres and normal[2] <= tres:
                   surf_vertical.append(k)
                   surf_vertical_tmp.append(k)
                   lcurve=cubit.get_relatives("surface",k,"curve")
                   list_curve_vertical=list_curve_vertical+list(lcurve)
                   command = "surface "+str(k)+" scheme submap"
                   cubit.cmd(command)
                   #cubit_error_stop(iproc,command,ner)
                   center_point = cubit.get_center_point("surface", k)
                elif normal[2] >= -1-tres and normal[2] <= -1+tres:
                   target_sweep=k
                   surf_or.append(k)
                   lcurve=cubit.get_relatives("surface",k,"curve")
                   list_curve_or=list_curve_or+list(lcurve)
                   command = "surface "+str(k)+" scheme pave"
                   cubit.cmd(command)
                   #cubit_error_stop(iproc,command,ner)
                   if id_vol == vol[0].ID:
                      top_bottom=k
                elif normal[2] >= 1-tres and normal[2] <= 1+tres:
                   if source_sweep == 0:
                      source_sweep = str(k)
                      source_sweep_add= ''
                   else:
                      source_sweep_add = str(k)
                   surf_or.append(k)
                   lcurve=cubit.get_relatives("surface",k,"curve")
                   list_curve_or=list_curve_or+list(lcurve)
                   if cfg.cut_outline and id_partitioner_surf in list_surface_transition and id_vol == vol[-1].ID :
                      if top_surface == 0:
                         top_surface=k
                      else:
                         top_surface_add=k
                         text='top add'+str(k)

                         area1=cubit.get_surface_area(top_surface)
                         area2=cubit.get_surface_area(top_surface_add)
                         if area1 > area2:
                            top_surface_add=top_surface
                            top_surface=k
                   elif id_vol == vol[-1].ID:
                      top_surface=k
                      top_surface_add=''
                   command = "surface "+str(k)+" scheme pave"
                   cubit.cmd(command)
                   #cubit_error_stop(iproc,command,ner)
                else:
                   text='normal = '+str(normal[2])
                   mpiprint(mpiflag,text)
                   text='treshold = '+str(tres)
                   mpiprint(mpiflag,text)
                   raise NameError, 'error assigning the curve' 
            for k in list_curve_or:
                length=cubit.get_curve_length(k)
                interval=int(2*round(.5*length/cfg.size,0))
                command = "curve "+str(k)+" interval "+str(interval)
                cubit.cmd(command)
                #cubit_error_stop(iproc,command,ner)
                command = "curve "+str(k)+" scheme equal"
                cubit.cmd(command)
                #cubit_error_stop(iproc,command,ner)
                icheck+=1
            for s in surf_vertical_tmp:
                lcurve=cubit.get_relatives("surface",s,"curve")
                interval_store=[]
                for k in lcurve:
                    interval_curve=cubit.get_mesh_intervals('curve',k)
                    if k not in list_curve_or:
                        command = "curve "+str(k)+" interval "+str(vol[l].intervalv)
                        cubit.cmd(command)
                        #cubit_error_stop(iproc,command,ner)
                        command = "curve "+str(k)+" scheme equal"
                        cubit.cmd(command)
                        #cubit_error_stop(iproc,command,ner)
                    else:
                        interval_store.append((k,interval_curve))
                if interval_store[0][1] == interval_store[1][1]:
                    pass
                else:
                    interval_min=min(interval_store[0][1],interval_store[1][1])
                    #I am sure that the minimum length is in the curve at the bottom so the min interval will be in vol[0] and the check will propagate correctly
                    command = "curve "+str(interval_store[0][0])+" "+str(interval_store[1][0])+" interval "+str(interval_min)
                    cubit.cmd(command)
                    #cubit_error_stop(iproc,command,ner)
                    command = "curve "+str(interval_store[0][0])+" "+str(interval_store[1][0])+" scheme equal"
                    cubit.cmd(command)
                    #cubit_error_stop(iproc,command,ner)
            command = "volume "+str(id_vol)+" scheme sweep source surface "+source_sweep+" "+source_sweep_add+" source target "+str(target_sweep)+ " rotate off"
            cubit.cmd(command)
            #cubit_error_stop(iproc,command,ner)
            ner=cubit.get_error_count()
            command = "volume "+str(id_vol)+" sweep smooth Auto"
            cubit.cmd(command)
            #cubit_error_stop(iproc,command,ner)
            #
        for id_volume in range(nvol-1,-1,-1):
            command = "mesh vol "+str(vol[id_volume].ID)
            cubit.cmd(command)
            #cubit_error_stop(iproc,command,ner)
        #
        #smoothing
        if  cfg.smoothing:
            cubitcommand= 'surf all smooth scheme laplacian '
            cubit.cmd(cubitcommand)
            #cubit_error_stop(iproc,cubitcommand,ner)
            cubitcommand= 'smooth surf all'
            cubit.cmd(cubitcommand)
            #cubit_error_stop(iproc,cubitcommand,ner)
            #
            cubitcommand= 'vol all smooth scheme laplacian '
            cubit.cmd(cubitcommand)
            #cubit_error_stop(iproc,cubitcommand,ner)
            cubitcommand= 'smooth vol all'
            cubit.cmd(cubitcommand)
            #cubit_error_stop(iproc,cubitcommand,ner)
        #
        #for nvol = 3 
        #
        #___________________________ interface 4
        #                 \*********
        #vol 2              \*******
        #_____________________\***** interface 3
        #
        #vol 1
        #___________________________ interface 2
        #
        #vol 0
        #___________________________ interface 1
        #
        #
        refinement(nvol,vol)
        if  cfg.coarsening_top_layer:
            from sets import set
            #
            cubitcommand= 'del mesh vol '+str(vol[-1].ID)+ ' propagate'
            cubit.cmd(cubitcommand)
            #cubit_error_stop(iproc,cubitcommand,ner)
            s1=set(list_curve_vertical)
            s2=set(cubit.parse_cubit_list("curve","in vol "+str(vol[-1].ID)))
            lc=list(s1 & s2).replace("["," ").replace("]"," ")
            cubitcommand= 'curve '+str(list(lc))+' interval '+str(cfg.actual_vertical_interval_top_layer)
            cubit.cmd(cubitcommand)
            #cubit_error_stop(iproc,cubitcommand,ner)
            cubitcommand= 'mesh vol '+str(vol[-1].ID)
            cubit.cmd(cubitcommand)
        #    
        #    
        if  cfg.refine_basin and id_partitioner_surf not in list_surface_nomesh:
            refinement_basin(iproc,id_partitioner_surf,list_surface_out,list_surface_transition,top_surface,top_surface_add)
        #
        #
        #    
        #save_view_file(iproc)
        #
        #
        n=cubit.get_sideset_id_list()
        if len(n) != 0:
            command = "del sideset all"
            cubit.cmd(command)
        n=cubit.get_block_id_list()
        if len(n) != 0:    
            command = "del block all"
            cubit.cmd(command)
        #
        #
        if cfg.refine_basin and cfg.interpolatehex and id_partitioner_surf in list_surface_transition:
                interpolate_basin(top_surface,top_surface_add)

        import boundary_definition
        entities=['face']
        if id_partitioner_surf in list_surface_transition:
            text='t'
        elif id_partitioner_surf in list_surface_out:
            text='out'
        elif id_partitioner_surf in list_surface_in:
            text='in'
        else:
            text=None
        boundary_definition.define_bc(entities,parallel=True,partitioning_type=text)

        #MPI scheme - checking boundaries
        if mpiflag and cfg.checkbound and len(list_surface_nomesh) != 0:
           import boundary_check
           boundary_check.bcheck_basin(surf_vertical,parallel_map)

        #
        #   
        if mpiflag and cfg.checkbound:
            if id_partitioner_surf not in list_surface_nomesh:
                command = "save as '"+cfg.output_dir+"/mesh_vol_ck_"+str(iproc)+".cub' overwrite"
                cubit.cmd(command)
                ##cubit_error_stop(iproc,command,ner)
                #cubit.cmd('export Genesis  "'+cfg.output_dir+'/mesh_vol_ck_'+str(iproc)+'.g" dimension 3 block all overwrite ')
                ##cubit_error_stop(iproc,command,ner)
            else:
                command = "save as '"+cfg.output_dir+"/NOmesh_vol_ck_"+str(iproc)+".cub' overwrite"
                cubit.cmd(command)
                ##cubit_error_stop(iproc,command,ner)
                #cubit.cmd('export Genesis  "'+cfg.output_dir+'/NOmesh_vol_ck_'+str(iproc)+'.g" dimension 3 block all overwrite ')
                ##cubit_error_stop(iproc,command,ner)
        else:
            print 'no mpi'
            command = "save as '"+cfg.output_dir+"/mesh_vol_nock"+str(iproc)+".cub' overwrite"
            cubit.cmd(command)
            #cubit.cmd('export Genesis  "'+cfg.output_dir+'mesh_vol_nock'+str(iproc)+'.g" dimension 3 block all overwrite ')
        #

        import quality_log
        max_skewness,min_length=quality_log.quality_log()
        count_hex=[cubit.get_hex_count()]
        count_node=[cubit.get_node_count()]
        max_skew=[(iproc,max_skewness)]
        min_l=[(iproc,min_length)]

        print 'waiting.....'
        mpi.barrier()
        total_min_l=mpi.gather(min_l)
        total_hex=mpi.gather(count_hex)        
        total_node=mpi.gather(count_node)      
        total_max_skew=mpi.gather(max_skew)    
        print iproc, ' ok '

        mpi.barrier()                          
        if iproc == 0:
            min_total_min_l=min([ms[1] for ms in total_min_l])
            max_total_max_skew=max([ms[1] for ms in total_max_skew])
            sum_total_node=sum(total_node)
            sum_total_hex=sum(total_hex)

            totstat_file=open(cfg.output_dir+'/totstat.log','w')
            text='hex total number,node total number,max skew, min length\n'
            totstat_file.write(text)
            print text
            text=str(sum_total_hex)+' , '+str(sum_total_node)+' , '+str(max_total_max_skew)+' , '+str(min_total_min_l)+'\n'
            totstat_file.write(text)
            print text
            totstat_file.write(str(total_max_skew))    
            totstat_file.close()