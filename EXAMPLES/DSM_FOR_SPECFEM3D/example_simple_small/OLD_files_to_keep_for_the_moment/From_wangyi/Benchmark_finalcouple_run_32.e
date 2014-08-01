+ /bin/bash -x /tmp/jobstart.12397
+ SCRIPT_PID=12411
+ set +x
++ id -u
+ '[' 24049 -ge 1000 ']'
+ MODULESHOME=/opt/Modules
+ export MODULESHOME
+ MODULEFILES=/opt/Modules/default/modulefiles
+ export MODULEFILES
+ '[' ccc:datadir/own:dfldatadir/own:bullxmpi/1.2.7.2:licsrv/intel:c/intel/14.0.3.174:c++/intel/14.0.3.174:fortran/intel/14.0.3.174:mkl/14.0.3.174:idb/14.0.3.174:intel/14.0.3.174 = '' ']'
+ ENV=/etc/profile.d/modules.sh
+ export ENV
+ BASH_ENV=/etc/profile.d/modules.sh
+ export BASH_ENV
+ FPATH=/opt/Modules/init/fpath
+ export FPATH
+ [[ hxB =~ i ]]
+ export module
+ [[ -s /opt/Modules/init/bash_completion ]]
+ [[ 4 -ge 3 ]]
+ [[ hxB =~ i ]]
+ '[' ccc:datadir/own:dfldatadir/own:bullxmpi/1.2.7.2:licsrv/intel:c/intel/14.0.3.174:c++/intel/14.0.3.174:fortran/intel/14.0.3.174:mkl/14.0.3.174:idb/14.0.3.174:intel/14.0.3.174 = '' ']'
+ set -x
+ cd /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small
+ declare -i NPROC NPROC_MINUS_ONE CPUS CHOICE MIDDLE
+ NPROC=32
+ CPUS=32
+ MPIRUN=ccc_mprun
+ OPTION=
+ NPROC_MINUS_ONE=32-1
++ pwd
+ flog_file=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/log.benchmark
+ PREFIX_MOVIE=velocity_Z_it
++ pwd
+ IN_MOVIE=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI/
++ pwd
+ OUT_MOVIE=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie/
+ /utils/DSM_FOR_SPECFEM3D/bin/xcreate_inputs_files
/tmp/jobstart.12397: line 106: /utils/DSM_FOR_SPECFEM3D/bin/xcreate_inputs_files: No such file or directory
+ echo '!!!!!!!!!!!!!!!!!! SHELLS STEP1 : fin de lecture parfile_for_benchmark !!!!!!!!!!!!!!!!'
+ CHOICE=3
+ '[' 3 -gt 3 ']'
+ '[' 3 -lt 1 ']'
+ MIDDLE=3
+ CHOICE=3
+ echo 'The value of CHOICE variable is' 3
+ echo 'The value of CHOICE variable is' 3
+ source params.in
++ HOME_SPECFEM3D=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d
++ BIN=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin
++ BINSEM=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin
++ SCRIPTS=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/shells
+++ pwd
++ DSM_tractions=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions/
++ OUT=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions/
++ REP=Tract/
+++ pwd
++ IN_DSM=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm
+++ pwd
++ MESH=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH
++ MODELE_1D=iasp91_dsm
+ source /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/shells/scrpits_specfem3D.sh
+ echo '!!!!!!!!!!!!!!!!!! SHELLS STEP3 : fin de lecture scrpits_specfem3D.sh !!!!!!!!!!!!!!!!'
+ '[' 3 -eq 1 ']'
+ '[' 3 -eq 2 ']'
+ source /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/shells/scripts_dsm_full.sh
+ echo
+ echo ' BENCHMARK RUN  '
+ echo
++ date
+ echo Fri Aug 1 16:19:28 CEST 2014
+ echo ''
++ date
+ echo Fri Aug 1 16:19:28 CEST 2014
+ run_create_specfem_databases
+ cp ParFileInterface bin/.
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xdecompose_mesh 32 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH OUTPUT_FILES/DATABASES_MPI/
+ mv Numglob2loc_elmn.txt /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/.
+ echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 1 !!!!!!!!!!!!!!!!'
+ pwd
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xgenerate_databases
+ echo ' create specfem3D data base'
++ date
+ echo Fri Aug 1 16:19:33 CEST 2014
+ exit 0
