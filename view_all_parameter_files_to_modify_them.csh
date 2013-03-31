#!/bin/csh

# this list of files is obtained using " find . -iname \*Par_file\* -exec ls -1 {} \; | grep -v Mesh_Par_file | grep -v Par_file_faults "
# and adding ../DATA/Par_file manually to the list

vi ../DATA/Par_file ./noise_tomography/DATA/Par_file_step3 ./noise_tomography/DATA/Par_file_step2 ./noise_tomography/DATA/Par_file_step1 ./tpv15/DATA/Par_file ./homogeneous_halfspace_HEX8/DATA/Par_file ./tpv16/DATA/Par_file ./tpv102/DATA/Par_file ./tpv5/DATA/Par_file ./splay_faults/DATA/Par_file ./homogeneous_poroelastic/DATA/Par_file ./tpv103/DATA/Par_file ./BENCHMARK_CLAERBOUT_ADJOINT/ACOUSTIC/DATA/Par_file ./BENCHMARK_CLAERBOUT_ADJOINT/ELASTIC/DATA/Par_file ./layered_halfspace/DATA/Par_file ./meshfem3D_examples/simple_model/DATA/Par_file ./meshfem3D_examples/socal1D/example_utm/Par_file_utm ./meshfem3D_examples/socal1D/DATA/Par_file ./meshfem3D_examples/many_interfaces/DATA/Par_file ./tomographic_model/DATA/Par_file ./homogeneous_halfspace_HEX27/DATA/Par_file ./waterlayered_halfspace/DATA/Par_file ./CPML/purely_elastic/HOMO8_FREE_SURFACE_with_CPML/DATA/Par_file ./CPML/purely_elastic/HOMO8_NOFREESURFACE_with_PML/DATA/Par_file ./CPML/purely_elastic/HOMO8_FREE_SURFACE_without_CPML/DATA/Par_file ./CPML/purely_acoustic/HOMO8_FREE_SURFACE/DATA/Par_file ./CPML/purely_acoustic/HOMO8_NOFREESURFACE/DATA/Par_file ./Mount_StHelens/DATA/Par_file

