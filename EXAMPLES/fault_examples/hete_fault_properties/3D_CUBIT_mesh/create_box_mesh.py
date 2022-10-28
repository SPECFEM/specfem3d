#!python
# Create 3D mesh files.
# Huihui Weng (Geoazur, 2018)

# ======================================================================
import numpy
import os
import sys
# Please set up the path for CUBIT (or Trelis) and GEOCUBIT in your system.
# It requires CUBIT version 14.0.5 or later
sys.path.append('/opt/linux64/Trelis/bin/')
sys.path.append('/opt/linux64/specfem3d/CUBIT_GEOCUBIT/')

import cubit
print "Init CUBIT..."
try:
    # print all the information to the screen.
     cubit.init([""])
    # stop all the outout information and warnings to the screen.
    #cubit.init(["-noecho","-nojournal","-information=off","-warning=off"])
except:
    pass
from geocubitlib import absorbing_boundary
from geocubitlib import save_fault_nodes_elements
from geocubitlib import cubit2specfem3d

#=====================================
#    Set up parameters             ===
#=====================================
# Please set up the mesh parametes in this section

# If DEBUG is True, then this script only create CUBIT script, otherwise create CUBIT script and mesh file.
# It is recommended to debug this script by GUI of CUBIT before to create Specfem3D mesh.
#DEBUG          = True
DEBUG          = False

# The dimension of model box (km)
Length         = 40
Width          = 30
Depth          = 20
# Move the box horizontally to be inside the range of interface and free surface (km)
Center_X       = 0
Center_Y       = 0

work_dir       = os.getcwd()
# If Interface is False, then use planar fault (given by the strike, dip, and dep). Otherwise run the scripts in ./Interface and give the path of the created interface (in the directory ./output)
# If Topography is False, then use planar surface. Otherwise run the scripts in ./Surface and give the path of the created planarsur (in the directory ./output)
Interface      =  False
Topography     =  False
Int_name       = work_dir + "/output/interface_sigma_1_inc_12.sat"
Top_name       = work_dir + "/output/surface_sigma_1_inc_12.sat"
Strike         = 90
Dip            = 90
# Indicating the vertical location of one reference point on fault, i.e., (0.0, 0.0, Dep)
Dep            = 0

# Uniform material properties.
vp  = 5770     # P wave speed (m/s)
vs  = 3330     # S wave speed (m/s)
rho = 2705     # density (g/m^3)
Q   = 13

# The mesh size (km). Smaller grid size can better sample curved geometries.
grid_size      = 1
# The mesh scheme: thex and map
# 1 -> Thex: firstly create a tetrahedral unstructured mesh, then convert into a hexahedral mesh (reduce the grid size by hal). This mesh scheme have good flexibility for curved geometries.
# 2 -> Map:  meshes all volumes with structured mesh of hexahedra. Before mesh by using this scheme, one needs to adjsut the domension of model box or move box horizontally to make all the surfaces have 4 sides. For example, if the fault cuts a volume that has a triangle surface, then the Map scheme doesn't work.
# Noted that the final mesh is hexahedral mesh
#mesh_scheme    = "thex"
mesh_scheme    = "map"
# The element type for hexahedral mesh: HEX8 or HEX27 (supported by Specfem3D)
# Higer order nodes can be moved to curved geometry by defaut, if set Node Constraint ON.
element_type = "HEX8"
#element_type = "HEX27"
# Refine the mesh of fault. fault_refine_numsplit=0 indicate  no refinement. fault_refine_numsplit (int) indicates how many times to subdivide the elements on fault.
# fault_refine_depth indicate the number of layers for refinement.
fault_refine_numsplit = 0
fault_refine_depth    = 5

# Set up the upper and lower depth of seimogenic zone. If Upper_cutoff>=0, then there is not upper seismogenic boundary and the fault cut through the free surface.
Upper_cutoff   = -3
#Upper_cutoff   =  0
Lower_cutoff   = -15

# The name of CUBIT script. One can run this script under the GUI of CUBIT for debuging. This python code will run this script without GUI.
journalFile    = "./output/Rupture_speed_animation.jou"
# The name (prefix name) of created mesh directory for Specfem3D. The full name is the prefix name + features of fault and free surface.
mesh_name      = "Rupture_speed_animation"
#=====================================


#==============================
#    Main code              ===
#==============================
# There is no need to change anything below. If you need to change something,
# please send me an email. I will try to make it more automatic.
#
print "Initial check..."
# Initial check
if(not os.path.isfile(Int_name) and Interface):
    print "The interface data does not exis!!! Please create it in ./Interface."
    exit()
elif(os.path.isfile(Int_name) and Interface):
    print "Using interface slab: ", Int_name
else:
    print "Using planar fault with strike: ", Strike, " dip: ", Dip, " depth(reference point): ", Dep

if(not os.path.isfile(Top_name) and Topography):
    print "The topography data does not exis!!! Please create it in ./Surface."
elif(os.path.isfile(Top_name) and Topography):
    print "Using topography: ", Top_name
else:
    print "Using planar topography."

# The name of output mesh file
if(Interface and Topography):
    output_mesh    = mesh_name + "_box_curvedfault_curvedtopo"
elif(not Interface and Topography):
    output_mesh    = mesh_name + "_box_planarfault" + "_strike_" + str(Strike) + "_dip_" + str(Dip) + "_depth_" + str(Dep) + "_curvedtopo"
elif(Interface and not Topography):
    output_mesh    = mesh_name + "_box_curvedfault_planarsur"
else:
    output_mesh    = mesh_name + "_box_planarfault" + "_strike_" + str(Strike) + "_dip_" + str(Dip) + "_depth_" + str(Dep) + "_planarsur"

# Add the info of upper boundary
if(Upper_cutoff>=0):
    output_mesh = output_mesh + "_cutsurf"
else:
    output_mesh = output_mesh + "_buried"

# Add the info of mesh scheme
output_mesh = output_mesh + "_size" + str(grid_size) + "_" + element_type

# Create the journal file for debuging
print "Create journal file..."
j = open(journalFile, 'w')
j.write("# Journal file formatting, etc.\n" + \
            "# ----------------------------------------------------------------------\n" + \
            "# Set units to SI.\n" + \
            "# ----------------------------------------------------------------------\n" \
            "${Units('si')}\n" + \
            "# Reset geometry.\n" + \
            "# ----------------------------------------------------------------------\n" \
            "reset\n")

j.write("# ----------------------------------------------------------------------\n" + \
            "# Create block\n" + \
            "# ----------------------------------------------------------------------\n")

j.write("${blockLength= %f *km}\n" % Length)
j.write("${blockWidth = %f *km}\n" % Width)
j.write("${blockHeight= %f *km}\n" % 2000)
j.write("brick x {blockLength} y {blockWidth} z {blockHeight}\n")
j.write("${idVol1=Id('volume')}\n")

j.write("# ----------------------------------------------------------------------\n" + \
            "# Move block\n" + \
            "# ----------------------------------------------------------------------\n")

j.write("${moveX= %f *km}\n" % Center_X)
j.write("${moveY= %f *km}\n" % Center_Y)
j.write("${moveZ= %f *km}\n" % 0)

j.write("volume {idVol1} move x {moveX} y {moveY} z {moveZ}\n")

if(Interface):
    j.write("# ----------------------------------------------------------------------\n" + \
            "# Import interface data.\n" + \
            "# ----------------------------------------------------------------------\n")
    j.write("import Acis '%s'\n" % Int_name)
    j.write("${idInt=Id('surface')}\n")
else:
    j.write("# ----------------------------------------------------------------------\n" + \
            "# Create planar interface.\n" + \
            "# ----------------------------------------------------------------------\n")
    j.write("create planar surface zplane\n")
    j.write("${idInt=Id('surface')}\n")
    j.write("rotate surface {idInt} about Y angle %f\n" % Dip)
    if(Strike != 0):
        j.write("rotate surface {idInt} about Z angle %f\n" % -Strike)
    j.write("surface {idInt} move z {%f}\n" % Dep)
if(Topography):
    j.write("# ----------------------------------------------------------------------\n" + \
            "# Import topography data\n" + \
            "# ----------------------------------------------------------------------\n")
    j.write("import Acis '%s'\n" % Top_name)
    j.write("${idSur=Id('surface')}\n")
else:
    j.write("# ----------------------------------------------------------------------\n" + \
            "# Create planar free surface.\n" + \
            "# ----------------------------------------------------------------------\n")
    j.write("create planar surface zplane\n")
    j.write("${idSur=Id('surface')}\n")

j.write("# ----------------------------------------------------------------------\n" + \
        "# Create bottom surface.\n" + \
        "# ----------------------------------------------------------------------\n")
j.write("create planar surface zplane offset {-%f *km}\n" % Depth)
j.write("${idBot=Id('surface')}\n")

j.write("# ----------------------------------------------------------------------\n" + \
        "# Webcut 1 block to 5 blocks.\n" + \
        "# ----------------------------------------------------------------------\n")
if(Upper_cutoff<0):
    j.write("webcut volume {idVol1} with sheet surface {idSur}\n")
    j.write("${idVol2=Id('volume')}\n")
    j.write("webcut volume {idVol2} with plane Zplane offset {%f *km}\n" % Upper_cutoff)
    j.write("${idVol3=Id('volume')}\n")
    j.write("webcut volume {idVol3} with plane Zplane offset {%f *km}\n" % Lower_cutoff)
    j.write("${idVol4=Id('volume')}\n")
    j.write("webcut volume {idVol4} with sheet surface {idBot}\n")
    j.write("${idVol5=Id('volume')}\n")
else:
    j.write("webcut volume {idVol1} with sheet surface {idSur}\n")
    j.write("${idVol3=Id('volume')}\n")
    j.write("webcut volume {idVol3} with plane Zplane offset {%f *km}\n" % Lower_cutoff)
    j.write("${idVol4=Id('volume')}\n")
    j.write("webcut volume {idVol4} with sheet surface {idBot}\n")
    j.write("${idVol5=Id('volume')}\n")


if(Interface):
    j.write("webcut volume {idVol3} with sheet surface {idInt}\n")
else:
    j.write("webcut volume {idVol3} with plane surface {idInt}\n")
j.write("${idVol6=Id('volume')}\n")
j.write("find surface overlap volume {idVol3} {idVol6}\n")
j.write("${idF1=GroupMemberId('surf_overlap','surface',0)}\n")
j.write("${idF2=GroupMemberId('surf_overlap','surface',1)}\n")
j.write("surface {idF1} name 'fault1'\n")

j.write("# ----------------------------------------------------------------------\n" + \
        "# Delete unused blocks and surfaces.\n" + \
        "# ----------------------------------------------------------------------\n")
j.write("delete surface all\n")
j.write("delete volume {idVol1} {idVol5} \n")

j.write("# ----------------------------------------------------------------------\n" + \
        "# imprint and merge.\n" + \
        "# ----------------------------------------------------------------------\n")
j.write("imprint all\n" + \
        "merge all\n")

j.write("# ----------------------------------------------------------------------\n" + \
            "# Generate the mesh.\n" + \
            "# ----------------------------------------------------------------------\n")
if(mesh_scheme == "thex"):
    j.write("volume all scheme TetMesh\n" + \
            "volume all size {%f*km}\n" % grid_size)
    j.write("mesh volume {idVol3} {idVol6}\n")
    j.write("mesh volume all\n")
    j.write("THex Volume all\n")
elif(mesh_scheme == "map"):
    j.write("volume all scheme map\n" + \
            "volume all size {%f*km}\n" % grid_size)
    j.write("mesh volume {idVol3} {idVol6}\n")
    j.write("mesh volume all\n")
else:
    print "Error mesh scheme!"
    exit()
if(fault_refine_numsplit > 0):
    j.write("refine surface fault1 NumSplit {0} depth {1}\n".format(fault_refine_numsplit,fault_refine_depth))

j.write("# ----------------------------------------------------------------------\n" + \
            "# Smooth mesh to improve quality.\n" + \
            "# ----------------------------------------------------------------------\n")
j.write("volume all smooth scheme condition number beta 2.0 cpu 4\n" + \
        "smooth volume all\n")

j.write("set unmerge Duplicate_mesh on\n")
j.write("unmerge surface fault1 only\n")
j.write("surface {idF2} name 'fault2'\n")
#j.write("unmerge curve in surface fault1\n")

j.write("# ----------------------------------------------------------------------\n" + \
            "# Seperate nodes on fault.\n" + \
            "# ----------------------------------------------------------------------\n")
j.write("set node constraint off\n")
j.write("node in surface fault1 move normal to surface fault1 distance {-0.01*m}\n")
j.write("node in surface fault2 move normal to surface fault2 distance {-0.01*m}\n")
j.write("compress all\n")
j.write("set node constraint on\n")

j.write("draw volume all\n")

j.write("# End of file\n")
j.close()

if(DEBUG):
   exit()

# ==================================================
#        Read the CUBIT journal and playback it.
# ==================================================
print "Playback journal file..."
with open(journalFile) as f:
    content = f.readlines()
for line in content:
    cubit.cmd(line)

# ==================================================
#         Save the mesh to txt files
#      This part is revised from the code of Specfem3D
# ==================================================
print ""
print "Convert mesh to Specfem-format..."
os.system('mkdir -p MESH')

## fault surfaces (up/down)
Au = [cubit.get_id_from_name("fault1")]
Ad = [cubit.get_id_from_name("fault2")]

#  FOR THE BULK (Seismic wave propagation part for SPECFEM3D)
entities=['face']
xmin,xmax,ymin,ymax,bottom,topo = absorbing_boundary.define_parallel_absorbing_surf()
# The above function is expected to obtain the correct topo and bottom surfaces, but for curved fault this function may fail (I don't know why).
# Here I try to obtain the topo and bottom surface automaticly in case the above function fails.
# If this part also doesn't work well (such as the topography has more than one surfaces), please setup bottom=[surface_list] and topo=[surface_list] manually.
list_surf    = cubit.parse_cubit_list("surface","all")
center_depth = []
for k in list_surf:
    center_depth.append(cubit.get_center_point("surface", k)[2])
bottom = [list_surf[numpy.argmin(center_depth)]]
topo   = [list_surf[numpy.argmax(center_depth)]]
#bottom = [surface_list]
#topo   = [surface_list]
if(len(bottom) == 0 or len(topo) == 0):
    print "Fail in obtaining the topo and bottom surfaces."
    print "Please change setup topo and bottom surfaces manually."
    exit()
print "Xmin surface list: ", xmin
print "Xmax surface list: ", xmax
print "Ymin surface list: ", ymin
print "Ymax surface list: ", ymax
print "Bott surface list: ", bottom
print "Topo surface list: ", topo

# define blocks
Vol_num = cubit.get_volume_count()
for i in range(Vol_num):
    cubit.cmd('block {0} hex in vol {0}'.format(i+1))
cubit.cmd('block 1003 face in surface ' + str(list(xmin)).replace("["," ").replace("]"," "))
cubit.cmd('block 1003 name "face_abs_xmin"')
cubit.cmd('block 1005 face in surface ' + str(list(xmax)).replace("["," ").replace("]"," "))
cubit.cmd('block 1005 name "face_abs_xmax"')
cubit.cmd('block 1004 face in surface ' + str(list(ymin)).replace("["," ").replace("]"," "))
cubit.cmd('block 1004 name "face_abs_ymin"')
cubit.cmd('block 1006 face in surface ' + str(list(ymax)).replace("["," ").replace("]"," "))
cubit.cmd('block 1006 name "face_abs_ymax"')
cubit.cmd('block 1002 face in surface ' + str(list(bottom)).replace("["," ").replace("]"," "))
cubit.cmd('block 1002 name "face_abs_bottom"')
cubit.cmd('block 1001 face in surface ' + str(list(topo  )).replace("["," ").replace("]"," "))
cubit.cmd('block 1001 name "face_abs_topo"')

# define boundaries as blocks
# Abandoned. This function creates blocks from 1 to n (first volume, then boundaries), but the function "export2SPECFEM3D" for creating hex27 need the boundaries to be defined from 1001 to 1006.
#absorbing_boundary.define_boundaries(entities, xmin,xmax,ymin,ymax,bottom,topo)

#### Define material properties for the 4 volumes ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')
for i in range(Vol_num):
    cubit.cmd('block {0}  name "elastic {0}" '.format(i+1))        # material region
    cubit.cmd('block {0} attribute count {1}'.format(i+1,6))
    cubit.cmd('block {0} attribute index 1 1'.format(i+1))
    cubit.cmd('block {0} attribute index 2 {1}'.format(i+1,vp))    # vp
    cubit.cmd('block {0} attribute index 3 {1}'.format(i+1,vs))    # vs
    cubit.cmd('block {0} attribute index 4 {1}'.format(i+1,rho))   # rho
    cubit.cmd('block {0} attribute index 5 {1}'.format(i+1,Q))     # Q flag (see constants.h: #IATTENUATION_ ... )
    cubit.cmd('block {0} attribute index 6 0'.format(i+1))         # q flag (see constants.h: iattenuation_ ... )

#### Export to SPECFEM3D format using cubit2specfem3d.py of GEOCUBIT
if(element_type == "HEX27"):
    cubit2specfem3d.export2SPECFEM3D('MESH',hex27=True)
else:
    cubit2specfem3d.export2SPECFEM3D('MESH')

# You need to create fault mesh file in the last, if using hex27.
faultA = save_fault_nodes_elements.fault_input(1,Au,Ad)

print "Save created mesh..."
# Save create directory as given name
os.system('rm -rf  output/' + output_mesh)
os.system('mv MESH output/' + output_mesh)

# End of script
