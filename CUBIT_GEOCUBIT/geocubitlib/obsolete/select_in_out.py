#############################################################################
# select_in_out.py                                                    
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

#####OBSOLETE
# use different scheme... numpy or cubit firetracing


def select_in_out(top_surf,curve_outline):
    import geocubitlib.menu as menu
    import geocubitlib.initializing as initializing
    #
    #
    mpiflag,iproc,numproc,mpi   = initializing.initializing_mpi()
    cubit                       = initializing.initializing_cubit()
    cfg                         = initializing.initializing_cfg()
    #
    store_node=[]
    command = "create planar surface with plane zplane offset 0"
    cubit.cmd(command)
    tmp_surface=cubit.get_last_id("surface")
    command = "project curve "+str(curve_outline)+" onto surf "+str(tmp_surface)
    cubit.cmd(command)
    outline_projected=cubit.get_last_id("curve")
    command = "group 'list_node' add node in surf "+str(top_surf)
    cubit.cmd(command)
    group=cubit.get_id_from_name("list_node")
    list_node=cubit.get_group_nodes(group)
    command = "delete group "+ str(group)
    cubit.cmd(command)
    #list_quad=cubit.get_surface_quads(top_surf)
    box=cubit.get_bounding_box("surface",tmp_surface)
    xmin=box[0]
    ymin=box[3]
    zmin=0
    for id_node in list_node:
        center_node=cubit.get_center_point("node",id_node)
        text='create curve location '+str(xmin)+' '+str(ymin)+' '+str(zmin)+' location '+str(center_node[0])+' '+str(center_node[1])+' 0'
        cubit.cmd(text)
        last_curve=cubit.get_last_id("curve")
        command = "project curve "+str(last_curve)+" on surf "+str(tmp_surface)
        cubit.cmd(command)
        command = "del curve "+str(last_curve)
        cubit.cmd(command)
        last_curve=cubit.get_last_id("curve")
        length=cubit.get_curve_length(last_curve)
        arc_vertex=cubit.get_last_id("vertex")
        command = "create vertex atintersection curve "+str(last_curve)+' '+str(outline_projected)
        cubit.cmd(command)
        last_vertex=cubit.get_last_id("vertex")
        icount=last_vertex-arc_vertex
        if icount != 0:
            for iv in range(arc_vertex+1,last_vertex+1):
                iv_coord=cubit.get_center_point("vertex",iv)
                distance=sqrt((xmin-iv_coord[0])**2+(ymin-iv_coord[1])**2)
                if distance > length:
                    icount=icount-1
            if icount%2 != 0:
                store_node.append(id_node)
        command = "delete curve "+str(last_curve)
        cubit.cmd(command)
        command = "del vertex all"
        cubit.cmd(command)
    command = "del surface "+str(tmp_surface)
    cubit.cmd(command)
    return store_node

def select_in_out_global(list_node,surf):
    from math import sqrt
    import geocubitlib.menu as menu
    import geocubitlib.initializing as initializing
    #
    #
    mpiflag,iproc,numproc,mpi   = initializing.initializing_mpi()
    cubit                       = initializing.initializing_cubit()
    cfg                         = initializing.initializing_cfg()
    #
    box=cubit.get_bounding_box("surface",surf)
    xmin=.5*(box[0]+box[1])
    ymin=.5*(box[3]+box[4])
    zmin=box[6]
    l=[]
    for id_node in list_node:
        center_node=cubit.get_center_point("node",id_node)
        text='create curve location '+str(xmin)+' '+str(ymin)+' '+str(zmin)+' location '+str(center_node[0])+' '+str(center_node[1])+' 0'
        cubit.cmd(text)
        last_curve=cubit.get_last_id("curve")
        last_v=cubit.get_last_id("vertex")
        last_s=cubit.get_last_id("surface")
        command='imprint volume in surface '+str(surf)+' with curve '+str(last_curve)
        cubit.cmd(command)
        last_s_after=cubit.get_last_id("surface")
        last_v_after=cubit.get_last_id("vertex")
        diff=last_v_after-last_v
        if diff % 2 != 0:
           l.append(id_node)
        command='delete curve '+str(last_curve)
        cubit.cmd(command)
        command='regularize vol in surf '+str(last_s) + ' to '+ str(last_s_after) 
        cubit.cmd(command)
    return l




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    