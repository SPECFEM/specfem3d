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

# The radius of the semi-sphere (km)
R_model        =  100
# The radius of the Cylinder that cut through both the free surface and fault (km)
R_cylinder     =  10

work_dir       = os.getcwd()
# If Interface is False, then use planar fault (given by the strike, dip, and dep). Otherwise run the scripts in ./Interface and give the path of the created interface (in the directory ./output)
# If Topography is False, then use planar surface. Otherwise run the scripts in ./Surface and give the path of the created planarsur (in the directory ./output)
Interface      = True
Topography     = True
Int_name       = work_dir + "/output/interface_sigma_1_inc_12.sat"
Top_name       = work_dir + "/output/surface_sigma_1_inc_12.sat"
Strike         = 230
Dip            = 70
# Indicating the vertical location of one reference point on fault, i.e., (0.0, 0.0, Dep)
Dep            = -5.7

# Uniform material properties.
vp  = 5770     # P wave speed (m/s)
vs  = 3330     # S wave speed (m/s)
rho = 2705     # density (g/m^3)
Q   = 13

# The mesh size (km). Smaller grid size can better sample curved geometries.
# fine_size is for the inner cylinder and coarse_size for the outside semisphere
fine_size        = 4
coarse_size      = 8
# The mesh scheme: thex
#  Thex: firstly create a tetrahedral unstructured mesh, then convert into a hexahedral mesh (reduce the grid size by hal). This mesh scheme have good flexibility for curved geometries.
#  Noted that the final mesh is hexahedral mesh
mesh_scheme    = "thex"
# The element type for hexahedral mesh: HEX8 or HEX27 (supported by Specfem3D)
# Higer order nodes can be moved to curved geometry by defaut, if set Node Constraint ON.
element_type = "HEX8"
#element_type = "HEX27"

# Set up the lower depth of seimogenic zone. The rupture can propogate to the free surface here.
Lower_cutoff   = -30
# The name of CUBIT script. One can run this script under the GUI of CUBIT for debuging. This python code will run this script without GUI.
journalFile    = "./output/Kumamoto.jou"
# The name (prefix name) of created mesh directory for Specfem3D. The full name is the prefix name + features of fault and free surface.
mesh_name      = "Kumamoto"
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
    output_mesh    = mesh_name + "_semisphere_curvedfault_curvedtopo"
elif(not Interface and Topography):
    output_mesh    = mesh_name + "_semisphere_planarfault" + "_strike_" + str(Strike) + "_dip_" + str(Dip) + "_depth_" + str(Dep) + "_curvedtopo"
elif(Interface and not Topography):
    output_mesh    = mesh_name + "_semisphere_curvedfault_planarsur"
else:
    output_mesh    = mesh_name + "_semisphere_planarfault" + "_strike_" + str(Strike) + "_dip_" + str(Dip) + "_depth_" + str(Dep) + "_planarsur"

# Add the info of mesh scheme
output_mesh = output_mesh + "_" + str(fine_size) + "_" + str(coarse_size) + "_" + element_type

# Vertical length of cylinder
L_cylinder = abs(2 * Lower_cutoff)

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
        "# Create a cylinder.\n" + \
        "# ----------------------------------------------------------------------\n")
j.write("create cylinder height {0} radius {1}\n".format("{"+str(L_cylinder)+"*km}",\
          "{"+str(R_cylinder)+"*km}"))
j.write("${idVol1=Id('volume')}\n")

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
    j.write("surface {idInt} move z {%f*km}\n" % Dep)
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
        "# Webcut blocks.\n" + \
        "# ----------------------------------------------------------------------\n")
j.write("webcut volume {idVol1} with sheet extended from surface {idSur}\n")
j.write("${idVol2=Id('volume')}\n")
j.write("webcut volume {idVol2} with sheet extended from surface {idInt}\n")
j.write("${idVol3=Id('volume')}\n")
j.write("# ----------------------------------------------------------------------\n" + \
        "# Find and name the fault surface.\n" + \
        "# ----------------------------------------------------------------------\n")
j.write("find surface overlap volume {idVol2} {idVol3}\n")
j.write("${idF1=GroupMemberId('surf_overlap','surface',0)}\n")
j.write("${idF2=GroupMemberId('surf_overlap','surface',1)}\n")
j.write("surface {idF1} name 'fault1'\n")

j.write("# ----------------------------------------------------------------------\n" + \
            "# Create semi-sphere\n" + \
            "# ----------------------------------------------------------------------\n")
j.write("create sphere radius {%f *km}\n" % R_model)
j.write("${idVol4=Id('volume')}\n")
j.write("webcut volume {idVol4} with sheet extended from surface {idSur}\n")
j.write("${idVol5=Id('volume')}\n")
j.write("${idround=Id('surface')}\n")
j.write("surface {idround} name 'spheresurf'\n")


j.write("# ----------------------------------------------------------------------\n" + \
            "# Substract the semi-spehere from the blocks that contain the fault\n" + \
            "# ----------------------------------------------------------------------\n")
j.write("subtract volume {idVol2} {idVol3} from volume {idVol5} keep\n")
j.write("${idVol6=Id('volume')}\n")

j.write("# ----------------------------------------------------------------------\n" + \
        "# Delete unused blocks and surfaces.\n" + \
        "# ----------------------------------------------------------------------\n")
j.write("delete surface all\n")
j.write("delete volume {idVol1} {idVol4} {idVol5} \n")

j.write("# ----------------------------------------------------------------------\n" + \
        "# imprint and merge.\n" + \
        "# ----------------------------------------------------------------------\n")
j.write("imprint all\n" + \
        "merge all\n")

j.write("# ----------------------------------------------------------------------\n" + \
            "# Generate the mesh.\n" + \
            "# ----------------------------------------------------------------------\n")
if(mesh_scheme == "thex"):
    j.write("volume all scheme TetMesh\n")
    j.write("volume {idVol2} {idVol3} size {%f*km}\n" % fine_size)
    j.write("mesh volume {idVol2} \n")
    j.write("mesh volume {idVol3}\n")
    j.write("volume {idVol6} size {%f*km}\n" % coarse_size)
    j.write("mesh volume {idVol6} \n")
    j.write("THex Volume all\n")
else:
    print "Error mesh scheme!"
    exit()

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

### Obtain the id of boundaries
# I define the original sphere surface as spheresurf. After webcut, CUBIT renames the new-cutted surface by adding @A, @B ...
SpheresurfID = [cubit.get_id_from_name("spheresurf@A")]
# Find the surface ID for the free surface
freesur_tolerance = 3e3
FreesurfID = []
list_surf=cubit.parse_cubit_list("surface","all")
for k in list_surf:
    center_point = cubit.get_center_point("surface", k)
    if abs(center_point[2]) <= freesur_tolerance:
         FreesurfID.append(k)

print SpheresurfID,FreesurfID
# define blocks
Vol_num = cubit.get_volume_count()
for i in range(Vol_num):
    cubit.cmd('block {0} hex in vol {0}'.format(i+1))
cubit.cmd('block 1000 face in surface ' + str(list(SpheresurfID)).replace("["," ").replace("]"," "))
cubit.cmd('block 1000 name "face_semisphere"')
cubit.cmd('block 1001 face in surface ' + str(list(FreesurfID)).replace("["," ").replace("]"," "))
cubit.cmd('block 1001 name "face_topo"')



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
