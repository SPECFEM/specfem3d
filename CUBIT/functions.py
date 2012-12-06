#!python

# Pg , functions Read_fault_nodes and nodes_coords_fault_open has been modified
#  to read properly fault_file_xxx.dat .
# Convention : Side up   = side 2   (+)
#            : Side down = side 1 (-)
import os

def read_nodes_coord_file(filein):
    file = open(filein,'r')
    nnodes = file.readline()
    nodxyzlist = file.readlines()
    x=[]
    y=[]
    z=[]
    node=[]
    for nodxyz in nodxyzlist:
        nodxyzsplit = nodxyz.split()
        node.append(int(nodxyzsplit[0]))
        x.append(float(nodxyzsplit[1]))
        y.append(float(nodxyzsplit[2]))
        z.append(float(nodxyzsplit[3]))
    file.close() 
    return node,x,y,z

def read_fault_nodes(filein):
    file = open(filein,'r')
    nodes_side1_side2 = file.readline()
    nside_du = nodes_side1_side2.split()
    nnodes_d = int(nside_du[0])
    nnodes_u = int(nside_du[1])
    nodes_u=[]
    nodes_d=[]
   
    for i in range(nnodes_d):
        hnodes = file.readline()
        hnodesplit = hnodes.split()
        nodes_d.append(int(hnodesplit[1]))
        nodes_d.append(int(hnodesplit[2]))
        nodes_d.append(int(hnodesplit[3]))
        nodes_d.append(int(hnodesplit[4]))

    for i in range(nnodes_u):
        hnodes = file.readline()
        hnodesplit = hnodes.split()
        nodes_u.append(int(hnodesplit[1]))
        nodes_u.append(int(hnodesplit[2]))
        nodes_u.append(int(hnodesplit[3]))
        nodes_u.append(int(hnodesplit[4]))

    nodes_d  = set(nodes_d)     # save unique fault nodes side u
    nodes_u  = set(nodes_u)     # save unique fault nodes side d
    ncommond = nodes_u&nodes_d  # glue nodes . 
    nodes_d  = nodes_d^ncommond # saving only split nodes side d
    nodes_u  = nodes_u^ncommond # saving only split nodes side u
    file.close()
    return nodes_u,nodes_d

# INPUTS :
#file_nodes_coord = 'MESH/nodes_coords_file'
#fault_file =       'MESH/fault_file_1.dat'
# OUTPUTS :
#name_out = 'MESH/nodes_coords_file_open_fault'
#fsideu   = 'MESH/fault_sideu.dat'
#fsided   = 'MESH/fault_sided.dat'
#nodes_coords_fault_open(file_nodes_coord,fault_file,name_out,fsideu,fsided,delta)
def nodes_coords_fault_open(file_nodes_coord,fault_file,name_out,fsideUP,fsideDW,delta): 
    x=[]
    y=[]
    z=[]
    node=[]
    nodes_u=[]
    nodes_d=[]
    txt = ''
    txt_fault = ''
    output_file = open(name_out,'w')
    fsideu = open(fsideUP,'w')
    fsided = open(fsideDW,'w')
    node,x,y,z = read_nodes_coord_file(file_nodes_coord)
    output_file.write('%10i\n' % len(node))
    node_u,node_d = read_fault_nodes(fault_file)
    
    for i in range(len(node)):
        if node[i] in node_u:
             z[i] = z[i] + delta/2.0
            # x[i] = x[i] + delta/2.0
             print "node_u",node[i]
             txt_fault=('%10i %20f %20f %20f\n') % (node[i],x[i],y[i],z[i])
             fsideu.write(txt_fault)
        if node[i] in node_d:
             z[i] = z[i] - delta/2.0
            # x[i] = x[i] - delta/2.0
             print "node_d",node[i]
             txt_fault=('%10i %20f %20f %20f\n') % (node[i],x[i],y[i],z[i])
             fsided.write(txt_fault)
        txt=('%10i %20f %20f %20f\n') % (node[i],x[i],y[i],z[i])
        output_file.write(txt)

    output_file.close()
    fsideu.close()
    fsided.close()
    return


delta = 0.5 # WARNING: Make sure that this delta is smaller than FAULT_GAP_TOLERANCE 
            # defined in decompose_mesh_SCOTH/fault_scotch.f90 (usually FAULT_GAP_TOLERANCE=1.0d0)
            # and larger than SMALLVAL_TOL defined in constants.h (usually SMALLVAL_TOL=1.d-10*dabs(UTM_X_MAX - UTM_X_MIN))

def m2km():
   name_in = 'MESH/nodes_coords_file'
   name_out = 'MESH/nodes_coords_file'
   txt =''
   km = 1000
   input_file  = open(name_in,'r')
   nnodes      = input_file.readline()
   nodes_list  = input_file.readlines()
   input_file.close()
   output_file = open(name_out,'w')
   nodes_num   = int(nnodes)
   output_file.write('%10i\n' % nodes_num)
   nodes       = []
   for node in nodes_list:
       nodes   = node.split()
       node_id = int(nodes[0])
       x= float(nodes[1])*km
       y= float(nodes[2])*km
       z= float(nodes[3])*km
       txt=('%10i %20f %20f %20f\n') % (node_id,x,y,z)
       output_file.write(txt)
   output_file.close()
   return
