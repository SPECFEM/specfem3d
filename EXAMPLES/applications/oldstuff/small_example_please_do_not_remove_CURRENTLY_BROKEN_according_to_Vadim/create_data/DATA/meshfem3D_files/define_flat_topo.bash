#
#        create topography free-surface for FK/Specfem hybrid method
#
#        z=0 is by convention at the bottom of layered model
#
#
#

dir=../../../../EXTERNAL_CODES_coupled_with_SPECFEM3D/AxiSEM_for_SPECFEM3D/utils_misc/


# number of point
nx=200
ny=200

# distance netween points
dx=1000.
dy=1000.

# first point
xmin=-1.
ymin=-1.

# elevation of topography
elevation=10000.

# flat topo
type_of_topo=fla   # gau or cos

echo $nx $ny > Par_topo
echo $dx $dy  >> Par_topo
echo $xmin $ymin  >> Par_topo
echo $type_of_topo >> Par_topo
echo $elevation >> Par_topo

$dir/create_interface.x<<EOF
Par_topo
EOF

awk '{print $3} ' interf.txt > interf.dat
