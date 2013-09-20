#############################################################################
# bc.py                                                    
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
def bcheck_regularmap():
    #
    import geocubitlib.initializing as initializing
    mpiflag,iproc,numproc,mpi   = initializing.initializing_mpi()
    from geocubitlib.mpi_geocubit import mpiprint
    #
    numpy                       = initializing.initializing_numpy()
    cubit                       = initializing.initializing_cubit()
    cfg                         = initializing.initializing_cfg()
    #
    if numproc != cfg.nproc_xi*cfg.nproc_eta:
        print 'check the number of processor, available '+str(numproc)+' cpus, requested '+str(cfg.nproc_xi*cfg.nproc_eta)+' cpus'
        import sys
        sys.exit()
    #
    from geocubitlib.utilities import geo2utm,cubit_error_stop
    ner=cubit.get_error_count()
    #
    from math import sqrt
    #
    ner=cubit.get_error_count()
    #
    x_slice=zeros([numproc],int)
    y_slice=zeros([numproc],int)
    for icpuy in range(0,cfg.nproc_eta): 
        for icpux in range (0,cfg.nproc_xi):
            iprocnum=icpuy*cfg.nproc_xi+icpux
            x_slice[iprocnum]=icpux
            y_slice[iprocnum]=icpuy
    #
    icpux=x_slice[iproc]
    icpuy=y_slice[iproc]
    #
    list_vol=cubit.parse_cubit_list("volume","all")
    nvol=len(list_vol)
    #
    #
    #                              
    zmax_box=cubit.get_total_bounding_box("volume",list_vol)[7]
    zmin_box=cubit.get_total_bounding_box("volume",list_vol)[6] #it is the z_min of the box ... box= xmin,xmax,d,ymin,ymax,d,zmin...    
    xmin_box=cubit.get_total_bounding_box("volume",list_vol)[0]
    xmax_box=cubit.get_total_bounding_box("volume",list_vol)[1]
    ymin_box=cubit.get_total_bounding_box("volume",list_vol)[3]
    ymax_box=cubit.get_total_bounding_box("volume",list_vol)[4]
    #
    command = "del sideset all"
    cubit.cmd(command)
    #
    cubitcommand= 'sideset '+str(3)+  ' surface with y_max <= '+str(ymin_box)
    cubit.cmd(cubitcommand)
    #
    cubitcommand= 'sideset '+str(4)+ ' surface with x_max <= '+str(xmin_box)
    cubit.cmd(cubitcommand)
    #
    cubitcommand= 'sideset '+str(5)+  ' surface with y_min >= '+str(ymax_box)
    cubit.cmd(cubitcommand)
    #
    cubitcommand= 'sideset '+str(6)+  ' surface with x_min >= '+str(xmax_box)
    cubit.cmd(cubitcommand)
    #
    cubit.cmd('group "xmin" add node in surface in sideset ' + str(4))
    cubit.cmd('group "xmax" add node in surface in sideset ' + str(6))
    cubit.cmd('group "ymin" add node in surface in sideset ' + str(3))
    cubit.cmd('group "ymax" add node in surface in sideset ' + str(5))
    # get the group ids
    group_xmin = cubit.get_id_from_name("xmin")
    group_xmax = cubit.get_id_from_name("xmax")
    group_ymin = cubit.get_id_from_name("ymin")
    group_ymax = cubit.get_id_from_name("ymax")
    # get the nodes in each group
    nodes_xmin = cubit.get_group_nodes(group_xmin)
    nodes_xmax = cubit.get_group_nodes(group_xmax)
    nodes_ymin = cubit.get_group_nodes(group_ymin)
    nodes_ymax = cubit.get_group_nodes(group_ymax)
    #
    cubit.cmd('group "ed" add edge in surf all')
    ie= cubit.get_id_from_name("ed")
    list_edge=cubit.get_group_edges(ie)
    len_edge=1.e9
    for edge_id in list_edge:
        len_edge_tmp=cubit.get_mesh_edge_length(edge_id)    
        if len_edge_tmp < len_edge: 
            len_edge = len_edge_tmp
    min_len=len_edge/2.
    max_len=len_edge/2.
    #
    i=0
    x_xmin=zeros([len(nodes_xmin)],float)
    y_xmin=zeros([len(nodes_xmin)],float)
    z_xmin=zeros([len(nodes_xmin)],float)
    for node_id in nodes_xmin:
        v = cubit.get_nodal_coordinates(node_id)
        x_xmin[i]=v[0]
        y_xmin[i]=v[1]
        z_xmin[i]=v[2]
        i+=1
    #
    #
    i=0
    x_xmax=zeros([len(nodes_xmax)],float)
    y_xmax=zeros([len(nodes_xmax)],float)
    z_xmax=zeros([len(nodes_xmax)],float)
    for node_id in nodes_xmax:
        v = cubit.get_nodal_coordinates(node_id)
        x_xmax[i]=v[0]
        y_xmax[i]=v[1]
        z_xmax[i]=v[2]
        i+=1
    #
    i=0
    x_ymin=zeros([len(nodes_ymin)],float)
    y_ymin=zeros([len(nodes_ymin)],float)
    z_ymin=zeros([len(nodes_ymin)],float)
    for node_id in nodes_ymin:
        v = cubit.get_nodal_coordinates(node_id)
        x_ymin[i]=v[0]
        y_ymin[i]=v[1]
        z_ymin[i]=v[2]
        i+=1
    #
    i=0
    x_ymax=zeros([len(nodes_ymax)],float)
    y_ymax=zeros([len(nodes_ymax)],float)
    z_ymax=zeros([len(nodes_ymax)],float)
    for node_id in nodes_ymax:
        v = cubit.get_nodal_coordinates(node_id)
        x_ymax[i]=v[0]
        y_ymax[i]=v[1]
        z_ymax[i]=v[2]
        i+=1
    #
    #
    #send rec the boundary nodes
    if icpux < cfg.nproc_xi-1:
        mpi.send(x_xmax,iproc+1,1)
        mpi.send(y_xmax,iproc+1,2)
        mpi.send(z_xmax,iproc+1,3)
    else:
        pass
    #
    if icpux > 0:
        x_xmin_ric,status=mpi.recv(iproc-1,1)
        y_xmin_ric,status=mpi.recv(iproc-1,2)
        z_xmin_ric,status=mpi.recv(iproc-1,3)
    else:
        x_xmin_ric=x_xmin
        y_xmin_ric=y_xmin
        z_xmin_ric=z_xmin
    #
    if icpuy < cfg.nproc_eta-1:
        mpi.send(x_ymax,iproc+cfg.nproc_xi,4)
        mpi.send(y_ymax,iproc+cfg.nproc_xi,5)
        mpi.send(z_ymax,iproc+cfg.nproc_xi,6)
    else:
        pass
    #
    if icpuy > 0:
        x_ymin_ric,status=mpi.recv(iproc-cfg.nproc_xi,4)
        y_ymin_ric,status=mpi.recv(iproc-cfg.nproc_xi,5)
        z_ymin_ric,status=mpi.recv(iproc-cfg.nproc_xi,6)
    else:
        x_ymin_ric=x_ymin
        y_ymin_ric=y_ymin
        z_ymin_ric=z_ymin
    #
    #
    #if icpux > 0:
    for node_id in nodes_xmin:
        n = cubit.get_nodal_coordinates(node_id)
        min_dist=1.e9 #change here
        for i in range(0,len(x_xmin_ric)):
            dx = n[0] - x_xmin_ric[i]
            dy = n[1] - y_xmin_ric[i]
            dz = n[2] - z_xmin_ric[i]
            dist = sqrt(dx*dx + dy*dy + dz*dz)
            if dist < min_dist:
               min_dist = dist
               arch_node = i
               if min_dist < cfg.precision: break
        if min_dist < min_len and min_dist > cfg.precision:
            command = "node " + str(node_id) + " move location position " + str(x_xmin_ric[arch_node]) +" "+ str(y_xmin_ric[arch_node]) + " "+str(z_xmin_ric[arch_node])
            cubit.cmd(command)
            command = "comment '"+str(list(n))+" -> "+command+"'"
            cubit.cmd(command)
            cubit_error_stop(iproc,command,ner)
        elif min_dist < cfg.precision and cfg.debug:
            command = "comment '"+"node "+str(node_id)+"below precision - no move"+"'"
            cubit.cmd(command)
        elif min_dist < cfg.precision:
            pass
        else:
            raise NameError, str(iproc)+' check boundaries failed, xmin'
    #
    #if icpuy > 0:
    for node_id in nodes_ymin:
        n = cubit.get_nodal_coordinates(node_id)
        min_dist=1.e9
        for i in range(0,len(x_ymin_ric)):
            dx = n[0] - x_ymin_ric[i]
            dy = n[1] - y_ymin_ric[i]
            dz = n[2] - z_ymin_ric[i]
            dist = sqrt(dx*dx + dy*dy + dz*dz)
            if dist < min_dist:
               min_dist = dist
               arch_node = i
               if min_dist < cfg.precision: break
        if min_dist < min_len and min_dist > cfg.precision:
            command = "node " + str(node_id) + " move location position " + str(x_ymin_ric[arch_node]) + " "+ str(y_ymin_ric[arch_node]) + " "+ str(z_ymin_ric[arch_node])
            cubit.cmd(command)
            cubit_error_stop(iproc,command,ner)
            command = "comment '"+str(list(n))+" -> "+command+"'"
            cubit.cmd(command)
        elif min_dist < cfg.precision and cfg.debug:
            command = "comment '"+"node "+str(node_id)+"below precision - no move"+"'"
            cubit.cmd(command)
        elif min_dist < cfg.precision:
            pass
        else:
            raise NameError, str(iproc)+' check boundaries failed, ymin'

def vert_surf_structure(iproc,surfv):
    import geocubitlib.initializing as initializing
    mpiflag,iproc,numproc,mpi   = initializing.initializing_mpi()
    from geocubitlib.mpi_geocubit import mpiprint
    #
    numpy                       = initializing.initializing_numpy()
    cubit                       = initializing.initializing_cubit()
    cfg                         = initializing.initializing_cfg()
    #
    #
    from math import sqrt
    from geocubitlib.utilities import geo2utm,cubit_error_stop
    #
    command = "comment '"+"vertical"+"'"
    cubit.cmd(command)
    command = "comment '"+str(surfv)+"'"
    cubit.cmd(command)
    #
    mapsurfvertical=[]
    mapsurfvertical.append(iproc)
    surf_properties=[]
    for s in surfv:
        p=cubit.get_center_point('surface',s)
        #cubit.cmd('del sideset all')
        command = 'group "n1" add node in surface '+str(s)
        cubit.cmd(command)
        group = cubit.get_id_from_name("n1")
        nodes = cubit.get_group_nodes(group)
        cubit.cmd('del group '+str(group))
        xnode=[]
        for node_id in nodes:
              v = cubit.get_nodal_coordinates(node_id)
              xnode.append([node_id,v])
        surf_properties=[s,p,xnode]
        mapsurfvertical.append(surf_properties)
    return mapsurfvertical
    #mapsurfvertical is: [ iproc number , [ ...[ surf number, [ xyz surface center point ] , [...[node number, xyz node]...]   ]...] ]


def cb(logfile,surf_local,surf_rec):
    import geocubitlib.initializing as initializing
    mpiflag,iproc,numproc,mpi   = initializing.initializing_mpi()
    from geocubitlib.mpi_geocubit import mpiprint
    #
    numpy                       = initializing.initializing_numpy()
    cubit                       = initializing.initializing_cubit()
    cfg                         = initializing.initializing_cfg()
    #
    #
    from math import sqrt
    from geocubitlib.utilities import geo2utm,cubit_error_stop
    #
    #
    #
    number_surf_l=len(surf_local)-1
    number_surf_r=len(surf_rec)-1
    #
    for s_l in surf_local[1:]:
        #
        flag=False
        for s_r in surf_rec[1:]:
            if abs(s_l[1][0]-s_r[1][0]) < cfg.precision and abs(s_l[1][1]-s_r[1][1]) < cfg.precision and abs(s_l[1][2]-s_r[1][2]) < cfg.precision:
               s=s_l[0]
               r=s_r[0]
               command = "comment ' surf "+str(s_l[0])+ "related to "+str(s_r[0])+"'"
               cubit.cmd(command)
               command = "comment '"+str(s_l[1])+"->"+str(s_r[1])+"'"
               cubit.cmd(command)
               #
               #
               if len(s_l[2]) != len(s_r[2]):
                    command = "comment 'ERROR: surf "+str(s_l[0])+ "has different number of nodes than surf "+str(s_r[0])+"'"
                    cubit.cmd(command)
                    command = "comment '"+str(len(s_l[2]))+"/"+str(len(s_r[2]))+"'"
                    cubit.cmd(command)
                    print>>logfile, '***********************************'
                    print>>logfile, 'surf local'+str(s_l[0])+' <-> '+'surf rec'+str(s_r[0])
                    print>>logfile, "comment 'ERROR: surf "+str(s_l[0])+ "has different number of nodes than surf "+str(s_r[0])+"'"
                    print>>logfile, "comment '"+str(len(s_l[2]))+"/"+str(len(s_r[2]))+"'"
                    print>>logfile, "comment '"+str(s_l[1])+"->"+str(s_r[1])+"'"
                    print>>logfile, 'local'
                    print>>logfile, s_l
                    print>>logfile, 'rec'
                    print>>logfile, s_r
               else:
                    p_l=s_l[2]
                    p_r=s_r[2]
                    flag=True
                    break
        if flag:
           cubit.cmd('group "ed" add edge in surf '+str(s_l[0]))
           ie= cubit.get_id_from_name("ed")
           list_edge=cubit.get_group_edges(ie)
           command = "del group "+str(ie)
           cubit.cmd(command)
           len_edge=1.e15
           for edge_id in list_edge:
               len_edge_tmp=cubit.get_mesh_edge_length(edge_id)    
               if len_edge_tmp < len_edge: 
                   len_edge = len_edge_tmp
           #
           min_len=len_edge/2.
           max_len=len_edge/2.
           #
           for point_l in p_l[1:]:
               flag_point=False
               min_dist=1.e9
               arch_node=point_l
               for point_r in p_r[1:]:
                   dx = point_l[1][0] - point_r[1][0]
                   dy = point_l[1][1] - point_r[1][1]
                   dz = point_l[1][2] - point_r[1][2]
                   dist = sqrt(dx*dx + dy*dy + dz*dz)
                   if dist < min_dist:
                      min_dist = dist
                      arch_node = point_r
                      if min_dist < cfg.precision: break
               point_r=arch_node
               if min_dist < min_len and min_dist > cfg.precision:
                  command = "node " + str(point_l[0]) + " move location position " + str(.5*(point_l[1][0] + point_r[1][0])) +" "+ str(.5*(point_l[1][1] + point_r[1][1])) + " "+ str(.5*(point_l[1][2] + point_r[1][2]))
                  cubit.cmd(command)
                  command = "comment '"+str(list(n))+" -> "+command+"'"
                  cubit.cmd(command)
               elif min_dist < cfg.precision and cfg.debug:
                  command = "comment '"+"node "+str(point_l[0])+" below precision "+str(min_dist)+" < "+str(cfg.precision)+" - no move"+"'"
                  cubit.cmd(command)
               else:
                  command = "comment '"+"node "+str(point_l[0])+" ERROR ["+str(point_l[1][0])+","+str(point_l[1][1])+","+str(point_l[1][2])+"] -> ["+str(point_r[1][0])+","+str(point_r[1][1])+","+str(point_r[1][2])+"] mindist="+str(min_dist)+" - ERROR NO MOVE"+"'"
                  cubit.cmd(command)





def bcheck_basin(surf_vertical,parallel_map):
    #
    import geocubitlib.initializing as initializing
    mpiflag,iproc,numproc,mpi   = initializing.initializing_mpi()
    from geocubitlib.mpi_geocubit import mpiprint
    #
    numpy                       = initializing.initializing_numpy()
    cubit                       = initializing.initializing_cubit()
    cfg                         = initializing.initializing_cfg()    
    #
    f=open(cfg.working+'mpi_tracking_'+str(iproc)+'.log','w')
    command = "comment 'check mpi boundaries proc: "+str(iproc)+"'"
    cubit.cmd(command)        
    #
    #
    surf_local=vert_surf_structure(iproc,surf_vertical)
    proc_adj=parallel_map[iproc].p_adj
    #
    mpi.barrier()
    #
    all_tmp=[]
    if iproc == 0:
        all_tmp.append(surf_local)
        for ip in range(1,numproc):
            command = "comment '"+"rec proc "+str(ip)+"'"
            cubit.cmd(command)
            s_tmp,status=mpi.recv(ip)
            all_tmp.append(s_tmp)
    else:
            mpi.send(surf_local,0)
    #all_surf_vertical=mpi.allgather(surf_local)
    
    mpi.barrier()
    all_surf_vertical=mpi.bcast(all_tmp)
    #
    #if iproc == 0:
    #    mpi.bcast(all_surf_vertical)
    #else:
    #    all_surf_vertical=mpi.bcast()
    #
    #print>>f,  all_surf_vertical
    #
    command = "comment '"+"checking..."+"'"
    cubit.cmd(command)
    for k in proc_adj:
        #for surf_rec in all_surf_vertical:
        #    if surf_rec[0] == k:
        surf_rec=all_surf_vertical[k]
        if surf_rec[0] != k:
           command = "comment '"+"ERROR: wrong vertical structure"+"'"
           cubit.cmd(command)
        else:
            command = "comment 'proc"+str(iproc)+"is checking with "+str(k)+"/"+str(proc_adj)+"'"
            cubit.cmd(command)
            cb(f,surf_local,surf_rec)

    