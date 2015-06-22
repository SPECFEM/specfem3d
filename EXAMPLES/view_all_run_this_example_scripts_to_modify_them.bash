#!/bin/bash

# this list of files is obtained using:    find . -iname run_this_example.sh -exec ls -1 {} \;

if [ -z "$EDITOR" ]
then
EDITOR=vi
fi

$EDITOR ./homogeneous_halfspace_HEX8_elastic_absorbing_Stacey_5sides/run_this_example.sh ./meshfem3D_examples/simple_model/run_this_example.sh ./meshfem3D_examples/many_interfaces/run_this_example.sh ./meshfem3D_examples/sep_bathymetry/run_this_example.sh ./meshfem3D_examples/socal1D/run_this_example.sh ./layered_halfspace/run_this_example.sh ./tomographic_model/run_this_example.sh ./waterlayered_halfspace/run_this_example.sh ./homogeneous_halfspace_HEX8_elastic_no_absorbing/run_this_example.sh ./Mount_StHelens/run_this_example.sh ./CPML_examples/homogeneous_halfspace_HEX8_acoustic_elastic_absorbing_CPML_5sides/run_this_example.sh ./CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_5sides/run_this_example.sh ./CPML_examples/homogeneous_halfspace_HEX8_acoustic_elastic_absorbing_CPML_6sides/run_this_example.sh ./CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_5sides/run_this_example.sh ./CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_6sides/run_this_example.sh ./CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_6sides/run_this_example.sh

