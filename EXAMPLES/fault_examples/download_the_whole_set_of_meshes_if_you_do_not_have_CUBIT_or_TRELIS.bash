#!/bin/bash

# this is optional, if you have CUBIT/TRELIS you can create the mesh files for all these fault_examples yourself by running the *.py files
# in CUBIT/TRELIS, otherwise you can download them from our external database

# this only needs to be done once (and for all)

echo " "
echo "downloading the whole set of mesh files for these fault_examples, which has a size of 680 MB; this may take a while..."
echo " "
wget http://data.geodynamics.org/specfem/specfem3d/all_mesh_files_for_fault_examples_splay_faults_tpv102_tpv103_tpv15_tpv16.tar.bz2

echo " "
echo "checking the checksum of th file downloaded..."
echo " (the result displayed should be 2ca335c7426432679fc0ddc040c49506, otherwise there was a transfer problem)"
echo " "
md5sum all_mesh_files_for_fault_examples_splay_faults_tpv102_tpv103_tpv15_tpv16.tar.bz2


echo " "
echo "uncompressing the file downloaded..."
echo " "
bunzip2 -f all_mesh_files_for_fault_examples_splay_faults_tpv102_tpv103_tpv15_tpv16.tar.bz2

