#!/bin/bash

# this list of files is obtained using " find . -iname \*Par_file\* -exec ls -1 {} \; | grep -v Mesh_Par_file | grep -v Par_file_faults | grep -v process_DATA_Par_files "
# and adding ../DATA/Par_file manually to the list

#
# Note that if you want to change the name of a parameter or anything else in all the Par_files below in an automatic way
# you can also use a "find" command combined with a "sed" command, here is an example:
#
#     find . -type f -name "*Par_file*" -exec sed -i 's/ATTENUATION_VISCOELASTIC_SOLID/ATTENUATION_VISCOELASTIC/g' {} \;
#
# However this will *not* make the change in ../DATA/Par_file, which you will need to do separately.
#

if [ -z "$EDITOR" ]
then
  EDITOR=vi
fi

$EDITOR ../DATA/Par_file ./homogeneous_poroelastic/DATA/Par_file ./layered_halfspace/DATA/Par_file ./meshfem3D_examples/socal1D/DATA/Par_file ./meshfem3D_examples/socal1D/example_utm/Par_file_utm ./meshfem3D_examples/many_interfaces/DATA/Par_file ./meshfem3D_examples/simple_model/DATA/Par_file ./meshfem3D_examples/sep_bathymetry/DATA/Par_file ./meshfem3D_examples/cavity/DATA/Par_file ./homogeneous_acoustic/DATA/Par_file ./homogeneous_halfspace_HEX8_elastic_absorbing_Stacey_5sides/DATA/Par_file ./BENCHMARK_CLAERBOUT_ADJOINT/ELASTIC/DATA/Par_file ./BENCHMARK_CLAERBOUT_ADJOINT/ACOUSTIC/DATA/Par_file ./CPML_examples/homogeneous_halfspace_HEX8_acoustic_elastic_absorbing_CPML_6sides/DATA/Par_file ./CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_6sides/DATA/Par_file ./CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_5sides/DATA/Par_file ./CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_5sides/DATA/Par_file ./CPML_examples/homogeneous_halfspace_HEX8_acoustic_elastic_absorbing_CPML_5sides/DATA/Par_file ./CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_6sides/DATA/Par_file ./homogeneous_halfspace_HEX27_elastic_no_absorbing/DATA/Par_file ./tomographic_model/DATA/Par_file ./homogeneous_halfspace_HEX8_elastic_no_absorbing/DATA/Par_file ./noise_tomography/DATA/Par_file_step2 ./noise_tomography/DATA/Par_file_step1 ./noise_tomography/DATA/Par_file_step3 ./Mount_StHelens/DATA/Par_file ./fault_examples/tpv103/DATA/Par_file ./fault_examples/tpv5/DATA/Par_file ./fault_examples/tpv16/DATA/Par_file ./fault_examples/splay_faults/DATA/Par_file ./fault_examples/tpv15/DATA/Par_file ./fault_examples/tpv102/DATA/Par_file ./waterlayered_halfspace/DATA/Par_file
