#!/bin/bash
#
# runs additional coverage examples
#

set -euo pipefail

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=$(pwd)
TESTID=${TESTID:-}
TESTCOV=${TESTCOV:-}

if [ "$TESTCOV" != "true" ]; then
  echo "TESTCOV=${TESTCOV} (not coverage), skipping run_coverage.sh"
  exit 0
fi

run_simple() {
  local rel_dir="$1"
  local nstep="$2"
  local model="${3:-}"

  echo "##################################################################"
  echo "${rel_dir}"
  echo

  cd "${WORKDIR}/${rel_dir}"

  # setup
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^NSTEP .*:NSTEP    = ${nstep}:" DATA/Par_file

  # noise
  if [ "${rel_dir}" == "EXAMPLES/applications/noise_tomography/" ]; then
    sed -i '10,$ d' NOISE_TOMOGRAPHY/S_squared    # truncates file, deletes all lines from line 10 till end
  fi
  # socal1D
  if [ "${rel_dir}" == "EXAMPLES/applications/meshfem3D_examples/socal1D/" ]; then
    if [ "$model" = "1d_socal" ]; then
      sed -i "s:^MODEL .*:MODEL    = 1d_socal:" DATA/Par_file
      rm -f REF_SEIS
      ln -s REF_SEIS.1d_socal REF_SEIS
    elif [ "$model" = "1d_prem" ]; then
      sed -i "s:^MODEL .*:MODEL    = 1d_prem:" DATA/Par_file
      rm -f REF_SEIS
      ln -s REF_SEIS.1d_prem REF_SEIS
    elif [ "$model" = "1d_cascadia" ]; then
      sed -i "s:^MODEL .*:MODEL    = 1d_cascadia:" DATA/Par_file
      rm -f REF_SEIS
      ln -s REF_SEIS.1d_cascadia REF_SEIS
    fi
  fi

  # run
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi

  # cleanup
  mv -v DATA/Par_file.org DATA/Par_file
  rm -rf OUTPUT_FILES/
  if [ -e DATABASES_MPI ]; then rm -rf DATABASES_MPI/; fi

  cd "$WORKDIR"
}

run_kernel() {
  local rel_dir="$1"
  local nstep="$2"

  echo "##################################################################"
  echo "${rel_dir} (kernel coverage)"
  echo

  cd "${WORKDIR}/${rel_dir}"

  # setup
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^NSTEP .*:NSTEP    = ${nstep}:" DATA/Par_file

  if [ "${rel_dir}" == "EXAMPLES/applications/homogeneous_acoustic/" ]; then
    sed -i "s:300:${nstep}:" run_this_example_kernel.sh
    sed -i "s:^t_start.*:t_start=-6.0:" create_adjoint_sources.sh
    sed -i "s:^t_end.*:t_end=-5.55:" create_adjoint_sources.sh
  fi

  # run
  ./run_this_example_kernel.sh
  if [[ $? -ne 0 ]]; then exit 1; fi

  # cleanup
  mv -v DATA/Par_file.org DATA/Par_file
  rm -rf OUTPUT_FILES/
  if [ -e DATABASES_MPI ]; then rm -rf DATABASES_MPI/; fi
  if [ -e SEM ]; then rm -rf SEM/; fi

  cd "$WORKDIR"
}


run_serial() {
  local rel_dir="$1"
  local nstep="$2"

  echo "##################################################################"
  echo "${rel_dir} (serial)"
  echo

  cd "${WORKDIR}/${rel_dir}"

  # setup
  cp -f DATA/Par_file DATA/Par_file.org
  sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = ${nstep}:" DATA/Par_file

  # meshfem setup
  if [ -e DATA/meshfem3D_files/Mesh_Par_file ]; then
    cp -f DATA/meshfem3D_files/Mesh_Par_file DATA/meshfem3D_files/Mesh_Par_file.org
    sed -i "s:^NPROC_XI .*:NPROC_XI    = 1:" DATA/meshfem3D_files/Mesh_Par_file
    sed -i "s:^NPROC_ETA .*:NPROC_ETA    = 1:" DATA/meshfem3D_files/Mesh_Par_file
  fi

  # run
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi

  # cleanup
  mv -v DATA/Par_file.org DATA/Par_file
  if [ -e DATA/meshfem3D_files/Mesh_Par_file.org ]; then mv -v DATA/meshfem3D_files/Mesh_Par_file.org DATA/meshfem3D_files/Mesh_Par_file; fi
  rm -rf OUTPUT_FILES/
  if [ -e DATABASES_MPI ]; then rm -rf DATABASES_MPI/; fi

  cd "$WORKDIR"
}

echo
echo "coverage run: TESTID=${TESTID}"
echo "work directory: ${WORKDIR}"
echo

# additional example tests (after base to avoid repeating code setup/configuration/compilation)
case "$TESTID" in
  0) # serial bunch
    run_serial "EXAMPLES/applications/homogeneous_halfspace/" 5
    run_serial "EXAMPLES/applications/meshfem3D_examples/simple_model/" 5
    run_serial "EXAMPLES/applications/Gmsh_simple_box_hex27/" 5
    ;;
  1) # parallel bunch 1
    run_simple "EXAMPLES/applications/homogeneous_halfspace_HEX27_elastic_no_absorbing/" 5
    run_simple "EXAMPLES/applications/homogeneous_poroelastic/" 5
    run_kernel "EXAMPLES/applications/homogeneous_acoustic/" 5
    ;;
  2) # parallel bunch 2
    run_simple "EXAMPLES/applications/CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_5sides/" 5
    run_simple "EXAMPLES/applications/CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_5sides/" 5
    run_simple "EXAMPLES/applications/noise_tomography/" 5
    run_simple "EXAMPLES/applications/tomographic_model/" 5
    run_simple "EXAMPLES/applications/layered_halfspace/" 5
    run_simple "EXAMPLES/applications/waterlayered_halfspace/" 5
    run_simple "EXAMPLES/applications/waterlayered_poroelastic/" 5
    ;;
  3) # parallel bunch 3
    run_simple "EXAMPLES/applications/fault_examples/tpv5/" 5
    run_simple "EXAMPLES/applications/small_example_coupling_FK_specfem/" 5
    run_simple "EXAMPLES/applications/Gmsh_simple_lddrk/" 5
    run_simple "EXAMPLES/applications/decompose_mesh_MPI/" 5
    run_simple "EXAMPLES/applications/small_adjoint_multiple_sources/" 5
    run_simple "EXAMPLES/applications/LTS_homogeneous_halfspace_HEX8/" 50
    ;;
  4) # parallel bunch 4
    run_simple "EXAMPLES/applications/meshfem3D_examples/socal1D/" 5 "1d_socal"
    run_simple "EXAMPLES/applications/meshfem3D_examples/socal1D/" 5 "1d_prem"
    run_simple "EXAMPLES/applications/meshfem3D_examples/socal1D/" 5 "1d_cascadia"
    run_simple "EXAMPLES/applications/meshfem3D_examples/cavity/" 5
    run_simple "EXAMPLES/applications/meshfem3D_examples/sep_bathymetry/" 5
    run_simple "EXAMPLES/applications/meshfem3D_examples/regular_element_mesh/" 5
    ;;
  5) # inversion
    echo "TESTID=4: no additional coverage examples (base inversion run already executed)"
    ;;
  6) # NGLL 6
    echo "TESTID=5: no additional coverage examples (base run already executed)"
    ;;
  *)
    echo "TESTID=${TESTID}: no additional coverage examples configured"
    ;;
esac

echo
echo "coverage examples done"
echo "$(date)"
echo
