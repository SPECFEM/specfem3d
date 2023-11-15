#!/bin/bash

# this list of files is obtained using:    find . -iname run_this_example.sh -exec ls -1 {} \;

if [ -z "$EDITOR" ]
then
EDITOR=vi
fi

# directories
dir=`pwd`
# changes to subdirectory EXAMPLES/ if called in root directory SPECFEM3D/
currentdir=`basename $dir`
echo "current directory: $currentdir"
if [ "$currentdir" == "SPECFEM3D" ]; then
cd EXAMPLES/
dir=`pwd`
else
echo "Please call this script from within the root directory SPECFEM3D/"; exit 1;
fi

$EDITOR ./applications/homogeneous_halfspace_HEX8_elastic_absorbing_Stacey_5sides/run_this_example.sh ./applications/meshfem3D_examples/simple_model/run_this_example.sh ./applications/meshfem3D_examples/many_interfaces/run_this_example.sh ./applications/meshfem3D_examples/sep_bathymetry/run_this_example.sh ./applications/meshfem3D_examples/socal1D/run_this_example.sh ./applications/layered_halfspace/run_this_example.sh ./applications/tomographic_model/run_this_example.sh ./applications/waterlayered_halfspace/run_this_example.sh ./applications/homogeneous_halfspace_HEX8_elastic_no_absorbing/run_this_example.sh ./applications/Mount_StHelens/run_this_example.sh ./applications/CPML_examples/homogeneous_halfspace_HEX8_acoustic_elastic_absorbing_CPML_5sides/run_this_example.sh ./applications/CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_5sides/run_this_example.sh ./applications/CPML_examples/homogeneous_halfspace_HEX8_acoustic_elastic_absorbing_CPML_6sides/run_this_example.sh ./applications/CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_5sides/run_this_example.sh ./applications/CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_6sides/run_this_example.sh ./applications/CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_6sides/run_this_example.sh

