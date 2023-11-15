#!/bin/bash

# this list of files is obtained using " find . -iname \*Par_file\* -exec ls -1 {} \; | /bin/grep -v Mesh_Par_file | /bin/grep -v Par_file_faults | /bin/grep -v process_DATA_Par_files | /bin/grep -v change_something_in_all_the_Par_files_automatically "

# ../DATA/Par_file is then automatically included in the list because it is a symbolic link to homogeneous_halfspace_HEX8_elastic_absorbing_Stacey_5sides/DATA/Par_file

#
# Note that if you want to change the name of a parameter or anything else in all the Par_files below in an automatic way
# you can also use a "find" command combined with a "sed" command, here is an example:
#
#     find . -type f -name "*Par_file*" -exec sed -i 's/ATTENUATION_VISCOELASTIC_SOLID/ATTENUATION_VISCOELASTIC/g' {} \;
#

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


$EDITOR ./applications/homogeneous_poroelastic/DATA/Par_file ./applications/homogeneous_poroelastic/DATA/Par_file_coarse ./applications/layered_halfspace/DATA/Par_file ./applications/meshfem3D_examples/socal1D/DATA/Par_file ./applications/meshfem3D_examples/socal1D/example_utm/Par_file_utm ./applications/meshfem3D_examples/many_interfaces/DATA/Par_file ./applications/meshfem3D_examples/simple_model/DATA/Par_file ./applications/meshfem3D_examples/sep_bathymetry/DATA/Par_file ./applications/meshfem3D_examples/cavity/DATA/Par_file ./applications/meshfem3D_examples/cavity/DATA/Par_file.96procs ./applications/homogeneous_acoustic/DATA/Par_file ./applications/small_example_coupling_axisem_specfem_script/Param_files/DATA/Par_file ./applications/small_example_coupling_axisem_specfem_script/Param_files/DATA/Par_file_one_proc ./applications/small_example_coupling_axisem_specfem_script/Param_files/DATA/Par_file_several_proc ./applications/small_example_coupling_axisem_specfem_script/Param_files_for_buried_box/DATA/Par_file ./applications/small_example_coupling_axisem_specfem_script/Param_files_for_buried_box/DATA/Par_file_one_proc ./applications/small_example_coupling_axisem_specfem_script/Param_files_for_buried_box/DATA/Par_file_several_proc ./applications/homogeneous_halfspace_HEX8_elastic_absorbing_Stacey_5sides/DATA/Par_file ./benchmarks/BENCHMARK_CLAERBOUT_ADJOINT/ELASTIC/DATA/Par_file ./benchmarks/BENCHMARK_CLAERBOUT_ADJOINT/ACOUSTIC/DATA/Par_file ./applications/CPML_examples/homogeneous_halfspace_HEX8_acoustic_elastic_absorbing_CPML_6sides/DATA/Par_file ./applications/CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_6sides/DATA/Par_file ./applications/CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_5sides/DATA/Par_file ./applications/CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_5sides/DATA/Par_file ./applications/CPML_examples/homogeneous_halfspace_HEX8_acoustic_elastic_absorbing_CPML_5sides/DATA/Par_file ./applications/CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_6sides/DATA/Par_file ./applications/Gmsh_simple_lddrk/DATA/Par_file ./applications/oldstuff/small_example_please_do_not_remove_CURRENTLY_BROKEN_according_to_Vadim/create_data/DATA/Par_file ./applications/oldstuff/small_example_please_do_not_remove_CURRENTLY_BROKEN_according_to_Vadim/fwi_fk/DATA/Par_file ./applications/homogeneous_halfspace_HEX27_elastic_no_absorbing/DATA/Par_file ./applications/tomographic_model/DATA/Par_file ./applications/homogeneous_halfspace_HEX8_elastic_no_absorbing/DATA/Par_file ./applications/small_example_coupling_axisem_specfem_matlab_gui_CURRENTLY_BROKEN_according_to_Vadim/Param_files/DATA/Par_file ./applications/small_example_coupling_axisem_specfem_matlab_gui_CURRENTLY_BROKEN_according_to_Vadim/Param_files/DATA/Par_file_one_proc ./applications/small_example_coupling_axisem_specfem_matlab_gui_CURRENTLY_BROKEN_according_to_Vadim/Param_files/DATA/Par_file_several_proc ./applications/noise_tomography/DATA/Par_file_step2 ./applications/noise_tomography/DATA/Par_file_step1 ./applications/noise_tomography/DATA/Par_file_step3 ./applications/Mount_StHelens/DATA/Par_file ./applications/small_example_coupling_FK_specfem/DATA/Par_file ./benchmarks/attenuation/viscoelastic/Par_file_attenuation_3D ./applications/fault_examples/tpv103/DATA/Par_file ./applications/fault_examples/tpv5/DATA/Par_file ./applications/fault_examples/tpv16/DATA/Par_file ./applications/fault_examples/splay_faults/DATA/Par_file ./applications/fault_examples/tpv15/DATA/Par_file ./applications/fault_examples/tpv102/DATA/Par_file ./applications/waterlayered_halfspace/DATA/Par_file
