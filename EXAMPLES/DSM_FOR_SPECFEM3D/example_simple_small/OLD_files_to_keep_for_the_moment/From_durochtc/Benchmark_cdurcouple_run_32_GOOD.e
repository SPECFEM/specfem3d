+ SCRIPT_PID=98173
+ /bin/bash -x /tmp/jobstart.98098
+ set +x
++ id -u
+ '[' 24229 -ge 1000 ']'
+ MODULESHOME=/opt/Modules
+ export MODULESHOME
+ MODULEFILES=/opt/Modules/default/modulefiles
+ export MODULEFILES
+ '[' ccc:datadir/own:dfldatadir/own:bullxmpi/1.2.7.2:licsrv/intel:c/intel/13.1.3.192:c++/intel/13.1.3.192:fortran/intel/13.1.3.192:mkl/13.1.3.192:idb/13.1.3.192:intel/13.1.3.192 = '' ']'
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
+ '[' ccc:datadir/own:dfldatadir/own:bullxmpi/1.2.7.2:licsrv/intel:c/intel/13.1.3.192:c++/intel/13.1.3.192:fortran/intel/13.1.3.192:mkl/13.1.3.192:idb/13.1.3.192:intel/13.1.3.192 = '' ']'
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
+ IN_MOVIE=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI
++ pwd
+ OUT_MOVIE=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xcreate_inputs_files
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
+ clean_and_make_dir
+ delete_directory_if_exist /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH
+ '[' '!' -d /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH ']'
+ delete_directory_if_exist OUTPUT_FILES
+ '[' '!' -d OUTPUT_FILES ']'
+ delete_directory_if_exist OUTPUT_FILES/DATABASES_MPI
+ '[' '!' -d OUTPUT_FILES/DATABASES_MPI ']'
+ delete_directory_if_exist OUTPUT_FILES/Tractions
+ '[' '!' -d OUTPUT_FILES/Tractions ']'
+ delete_directory_if_exist bin
+ '[' '!' -d bin ']'
+ mkdir -p /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH
+ mkdir -p OUTPUT_FILES/
+ mkdir -p OUTPUT_FILES/DATABASES_MPI/
+ mkdir -p OUTPUT_FILES/Tractions/
+ mkdir bin/
mkdir: cannot create directory `bin/': File exists
+ clean_and_make_dir_dsm
+ delete_directory_if_exist /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions/
+ '[' '!' -d /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions/ ']'
+ delete_directory_if_exist Tract/
+ '[' '!' -d Tract/ ']'
+ mkdir -p /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions/
+ mkdir -p Tract/
+ mv input_dsm_for_write_coef /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_for_coef
+ mv input_dsm_for_read_xmin /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_stxmin
+ mv input_dsm_for_read_xmax /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_stxmax
+ mv input_dsm_for_read_ymin /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_stymin
+ mv input_dsm_for_read_ymax /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_stymax
+ mv input_dsm_for_read_zmin /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_stzmin
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/iasp91 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/.
+ echo
+ echo ' BENCHMARK RUN  '
+ echo
++ date
+ echo Thu Jul 31 16:21:55 CEST 2014
+ echo ''
++ date
+ echo Thu Jul 31 16:21:55 CEST 2014
+ run_create_specfem_databases
+ cp ParFileInterface bin/.
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xdecompose_mesh 32 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH OUTPUT_FILES/DATABASES_MPI/
+ mv Numglob2loc_elmn.txt /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/.
+ echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 1 !!!!!!!!!!!!!!!!'
+ pwd
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xgenerate_databases
+ echo ' create specfem3D data base'
++ date
+ echo Thu Jul 31 16:22:00 CEST 2014
+ echo ''
+ echo ' create traction database'
++ date
+ echo Thu Jul 31 16:22:00 CEST 2014
+ run_create_tractions_for_specfem
+ cp ParFileInterface bin/.
+ echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 2 !!!!!!!!!!!!!!!!'
+ pwd
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xread_absorbing_interfaces
+ cp out_read.txt bin/
++ date
+ echo Thu Jul 31 16:24:23 CEST 2014
+ echo ''
+ echo ' simulation'
++ date
+ echo Thu Jul 31 16:24:23 CEST 2014
+ run_simu
+ echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 3 !!!!!!!!!!!!!!!!'
+ pwd
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xspecfem3D
+ cp out_read.txt bin/
+ echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 4 !!!!!!!!!!!!!!!!'
+ pwd
++ date
+ echo Thu Jul 31 16:27:16 CEST 2014
+ exit 0
