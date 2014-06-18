#!/bin/csh

# this list of files is obtained using " find . -iname \*proces\* -exec ls -1 {} \;

vi \
./layered_halfspace/process.sh \
./homogeneous_halfspace_HEX8_elastic_absorbing_Stacey_5sides/process.sh \
./BENCHMARK_CLAERBOUT_ADJOINT/ACOUSTIC/postprocessing.f90 \
./BENCHMARK_CLAERBOUT_ADJOINT/ELASTIC/postprocessing.f90 \
./waterlayered_halfspace/process.sh \
./Mount_StHelens/process.sh \
./noise_tomography/pre-processing.sh \
./CPML_examples/homogeneous_halfspace_HEX8_acoustic_elastic_absorbing_CPML_6sides/process.sh \
./CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_5sides/process.sh \
./CPML_examples/homogeneous_halfspace_HEX8_acoustic_elastic_absorbing_CPML_5sides/process.sh \
./CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_6sides/process.sh \
./CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_6sides/process.sh \
./CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_5sides/process.sh \
./meshfem3D_examples/many_interfaces/go_process_pbs.bash \
./meshfem3D_examples/many_interfaces/process.sh \
./meshfem3D_examples/simple_model/process.sh \
./meshfem3D_examples/socal1D/process.sh \
./homogeneous_halfspace_HEX8_elastic_no_absorbing/process.sh \
./tomographic_model/process.sh
