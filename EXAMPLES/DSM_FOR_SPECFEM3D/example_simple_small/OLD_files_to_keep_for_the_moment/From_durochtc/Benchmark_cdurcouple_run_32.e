+ SCRIPT_PID=78383
+ set +x
+ /bin/bash -x /tmp/jobstart.78369
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
+ source /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/shells/scripts_specfem3D.sh
+ echo '!!!!!!!!!!!!!!!!!! SHELLS STEP3 : fin de lecture scripts_specfem3D.sh !!!!!!!!!!!!!!!!'
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
+ mkdir -p OUTPUT_FILES
+ mkdir -p OUTPUT_FILES/DATABASES_MPI
+ mkdir -p OUTPUT_FILES/Tractions
+ mkdir bin
mkdir: cannot create directory `bin': File exists
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
+ echo Wed Jul 30 22:10:24 CEST 2014
+ run_create_mesh
++ pwd
+ current_dir=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small
+ cp ParFileMeshChunk /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/iasp91_dsm /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/.
+ cd /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xmesh_chunk_vm
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/model_1D.in ../DATA/.
+ cd /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small
+ cd /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions/
+ make_dir_exp
+ mkdir log
mkdir: cannot create directory `log': File exists
+ mkdir Displacement
mkdir: cannot create directory `Displacement': File exists
+ mkdir Stress
mkdir: cannot create directory `Stress': File exists
+ mkdir ascii
mkdir: cannot create directory `ascii': File exists
+ copy_input_files_exp
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/FrqsMpi.txt /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_for_coef /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//inputIASP.infTra
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/iasp91_dsm /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/st /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/recdepth /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//.
+ compute_exp_coeff
+ echo 32
+ echo
+ echo /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_write_ceof_mpi_SH
+ echo /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_write_ceof_mpi_PSV
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_write_ceof_mpi_SH
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_write_ceof_mpi_PSV
+ run_dsm_traction
+ cd /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions/
+ make_dir_faces
+ mkdir -p STXMIN/Displacement
+ mkdir -p STXMIN/Stress
+ mkdir -p STXMIN/log
+ mkdir -p STXMIN/out
+ mkdir -p STXMAX/Displacement
+ mkdir -p STXMAX/Stress
+ mkdir -p STXMAX/log
+ mkdir -p STXMAX/out
+ mkdir -p STYMIN/Displacement
+ mkdir -p STYMIN/Stress
+ mkdir -p STYMIN/log
+ mkdir -p STYMIN/out
+ mkdir -p STYMAX/Displacement
+ mkdir -p STYMAX/Stress
+ mkdir -p STYMAX/log
+ mkdir -p STYMAX/out
+ mkdir -p STZMIN/Displacement
+ mkdir -p STZMIN/Stress
+ mkdir -p STZMIN/log
+ mkdir -p STZMIN/out
+ copy_input_files_faces
+ STD=STXMIN
+ stf=stxmin
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/FrqsMpi.txt STXMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/Double_para.txt STXMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/iasp91_dsm STXMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_stxmin STXMIN/inputIASP.infTra
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/stxmin STXMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/recdepth STXMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/OrigRepSpecfm STXMIN/.
+ STD=STXMAX
+ stf=stxmax
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/FrqsMpi.txt STXMAX/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/Double_para.txt STXMAX/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/iasp91_dsm STXMAX/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_stxmax STXMAX/inputIASP.infTra
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/stxmax STXMAX/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/recdepth STXMAX/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/OrigRepSpecfm STXMAX/.
+ STD=STYMIN
+ stf=stymin
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/FrqsMpi.txt STYMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/Double_para.txt STYMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/iasp91_dsm STYMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_stymin STYMIN/inputIASP.infTra
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/stymin STYMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/recdepth STYMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/OrigRepSpecfm STYMIN/.
+ STD=STYMAX
+ stf=stymax
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/FrqsMpi.txt STYMAX/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/Double_para.txt STYMAX/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/iasp91_dsm STYMAX/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_stymax STYMAX/inputIASP.infTra
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/stymax STYMAX/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/recdepth STYMAX/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/OrigRepSpecfm STYMAX/.
+ STD=STZMIN
+ stf=stzmin
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/FrqsMpi.txt STZMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/Double_para.txt STZMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/iasp91_dsm STZMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/input_dsm/inputIASP.infTra_stzmin STZMIN/inputIASP.infTra
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/stzmin STZMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/recdepth STZMIN/.
+ cp /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/OrigRepSpecfm STZMIN/.
+ echo
+ echo ' FACE xmin'
++ date
+ echo Wed Jul 30 22:12:23 CEST 2014
+ cd STXMIN
+ read_exp_vert inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_MPI_read_ceof_vert_SH inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_MPI_read_ceof_vert_PSV inputIASP.infTra
++ date
+ echo Wed Jul 30 22:25:51 CEST 2014
+ fft_vert inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/fft_face_vert_full inputIASP.infTra
++ date
+ echo Wed Jul 30 22:26:06 CEST 2014
+ change_format_vertical
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/ChangeFormat
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/ChangeFormat_disp
++ date
+ echo Wed Jul 30 22:26:44 CEST 2014
+ echo
+ echo 'FACE xmax'
++ date
+ echo Wed Jul 30 22:26:44 CEST 2014
+ cd ../STXMAX
+ read_exp_vert inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_MPI_read_ceof_vert_SH inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_MPI_read_ceof_vert_PSV inputIASP.infTra
++ date
+ echo Wed Jul 30 22:33:38 CEST 2014
+ fft_vert inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/fft_face_vert_full inputIASP.infTra
++ date
+ echo Wed Jul 30 22:34:26 CEST 2014
+ change_format_vertical
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/ChangeFormat
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/ChangeFormat_disp
++ date
+ echo Wed Jul 30 22:35:25 CEST 2014
+ echo
+ echo 'FACE ymin'
++ date
+ echo Wed Jul 30 22:35:26 CEST 2014
+ cd ../STYMIN
+ read_exp_vert inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_MPI_read_ceof_vert_SH inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_MPI_read_ceof_vert_PSV inputIASP.infTra
++ date
+ echo Wed Jul 30 22:46:32 CEST 2014
+ fft_vert inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/fft_face_vert_full inputIASP.infTra
++ date
+ echo Wed Jul 30 22:47:02 CEST 2014
+ change_format_vertical
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/ChangeFormat
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/ChangeFormat_disp
++ date
+ echo Wed Jul 30 22:48:09 CEST 2014
+ echo
+ echo 'FACE ymax'
++ date
+ echo Wed Jul 30 22:48:09 CEST 2014
+ cd ../STYMAX
+ read_exp_vert inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_MPI_read_ceof_vert_SH inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_MPI_read_ceof_vert_PSV inputIASP.infTra
++ date
+ echo Wed Jul 30 22:53:43 CEST 2014
+ fft_vert inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/fft_face_vert_full inputIASP.infTra
++ date
+ echo Wed Jul 30 22:54:16 CEST 2014
+ change_format_vertical
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/ChangeFormat
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/ChangeFormat_disp
++ date
+ echo Wed Jul 30 22:55:23 CEST 2014
+ echo
+ echo 'FACE zmin'
++ date
+ echo Wed Jul 30 22:55:23 CEST 2014
+ cd ../STZMIN
+ read_exp_zmin inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_MPI_read_zmin_SH inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xTraPSV_MPI_read_zmin_PSV inputIASP.infTra
++ date
+ echo Wed Jul 30 23:00:01 CEST 2014
+ fft_zmin inputIASP.infTra
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/TraFFT_MPI_face_zmin_full inputIASP.infTra
++ date
+ echo Wed Jul 30 23:01:14 CEST 2014
+ change_format_zmin
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/ChangeFormat_zmin
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/ChangeFormat_zmin_disp
++ date
+ echo Wed Jul 30 23:02:19 CEST 2014
+ cd ../..
+ move_output
+ mkdir -p Tract/
+ mv /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//STXMIN/Trac.bin Tract//tractxmin.bin
+ mv /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//STXMIN/Disp.bin Tract//velxmin.bin
+ mv /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//STXMAX/Trac.bin Tract//tractxmax.bin
+ mv /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//STXMAX/Disp.bin Tract//velxmax.bin
+ mv /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//STYMIN/Trac.bin Tract//tractymin.bin
+ mv /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//STYMIN/Disp.bin Tract//velymin.bin
+ mv /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//STYMAX/Trac.bin Tract//tractymax.bin
+ mv /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//STYMAX/Disp.bin Tract//velymax.bin
+ mv /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//STZMIN/Trac.bin Tract//tractzmin.bin
+ mv /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/DSM_tractions//STZMIN/Disp.bin Tract//velzmin.bin
+ echo ''
++ date
+ echo Wed Jul 30 23:02:19 CEST 2014
+ run_create_specfem_databases
+ cp ParFileInterface bin/.
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xdecompose_mesh 32 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH OUTPUT_FILES/DATABASES_MPI/
+ mv Numglob2loc_elmn.txt /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/MESH/.
+ echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 1 !!!!!!!!!!!!!!!!'
+ pwd
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xgenerate_databases
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
error opening database proc######_external_mesh.bin
+ echo ' create specfem3D data base'
++ date
+ echo Wed Jul 30 23:02:25 CEST 2014
+ echo ''
+ echo ' create traction database'
++ date
+ echo Wed Jul 30 23:02:25 CEST 2014
+ run_create_tractions_for_specfem
+ cp ParFileInterface bin/.
+ echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 2 !!!!!!!!!!!!!!!!'
+ pwd
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/bin/xread_absorbing_interfaces
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B8708D36A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B8708D354B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B87080F519E  Unknown               Unknown  Unknown
libifcore.so.5     00002B8708064C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B870806414D  Unknown               Unknown  Unknown
libifcore.so.5     00002B870809A340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002ADB302DAA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002ADB302D94B6  Unknown               Unknown  Unknown
libifcore.so.5     00002ADB2F69919E  Unknown               Unknown  Unknown
libifcore.so.5     00002ADB2F608C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002ADB2F60814D  Unknown               Unknown  Unknown
libifcore.so.5     00002ADB2F63E340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002AF67BE30A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002AF67BE2F4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002AF67B1EF19E  Unknown               Unknown  Unknown
libifcore.so.5     00002AF67B15EC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002AF67B15E14D  Unknown               Unknown  Unknown
libifcore.so.5     00002AF67B194340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002AEAA86A3A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002AEAA86A24B6  Unknown               Unknown  Unknown
libifcore.so.5     00002AEAA7A6219E  Unknown               Unknown  Unknown
libifcore.so.5     00002AEAA79D1C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002AEAA79D114D  Unknown               Unknown  Unknown
libifcore.so.5     00002AEAA7A07340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002ADD927E3A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002ADD927E24B6  Unknown               Unknown  Unknown
libifcore.so.5     00002ADD91BA219E  Unknown               Unknown  Unknown
libifcore.so.5     00002ADD91B11C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002ADD91B1114D  Unknown               Unknown  Unknown
libifcore.so.5     00002ADD91B47340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002AD7B2900A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002AD7B28FF4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002AD7B1CBF19E  Unknown               Unknown  Unknown
libifcore.so.5     00002AD7B1C2EC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002AD7B1C2E14D  Unknown               Unknown  Unknown
libifcore.so.5     00002AD7B1C64340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B4256F0EA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B4256F0D4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B42562CD19E  Unknown               Unknown  Unknown
libifcore.so.5     00002B425623CC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B425623C14D  Unknown               Unknown  Unknown
libifcore.so.5     00002B4256272340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002BA9AF262A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002BA9AF2614B6  Unknown               Unknown  Unknown
libifcore.so.5     00002BA9AE62119E  Unknown               Unknown  Unknown
libifcore.so.5     00002BA9AE590C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002BA9AE59014D  Unknown               Unknown  Unknown
libifcore.so.5     00002BA9AE5C6340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B67857E2A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B67857E14B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B6784BA119E  Unknown               Unknown  Unknown
libifcore.so.5     00002B6784B10C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B6784B1014D  Unknown               Unknown  Unknown
libifcore.so.5     00002B6784B46340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002AC11E8A7A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002AC11E8A64B6  Unknown               Unknown  Unknown
libifcore.so.5     00002AC11DC6619E  Unknown               Unknown  Unknown
libifcore.so.5     00002AC11DBD5C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002AC11DBD514D  Unknown               Unknown  Unknown
libifcore.so.5     00002AC11DC0B340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B69B0A5BA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B69B0A5A4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B69AFE1A19E  Unknown               Unknown  Unknown
libifcore.so.5     00002B69AFD89C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B69AFD8914D  Unknown               Unknown  Unknown
libifcore.so.5     00002B69AFDBF340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002AE54890BA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002AE54890A4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002AE547CCA19E  Unknown               Unknown  Unknown
libifcore.so.5     00002AE547C39C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002AE547C3914D  Unknown               Unknown  Unknown
libifcore.so.5     00002AE547C6F340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B551AB75A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B551AB744B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B5519F3419E  Unknown               Unknown  Unknown
libifcore.so.5     00002B5519EA3C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B5519EA314D  Unknown               Unknown  Unknown
libifcore.so.5     00002B5519ED9340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B67BD578A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B67BD5774B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B67BC93719E  Unknown               Unknown  Unknown
libifcore.so.5     00002B67BC8A6C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B67BC8A614D  Unknown               Unknown  Unknown
libifcore.so.5     00002B67BC8DC340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002ABE66600A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002ABE665FF4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002ABE659BF19E  Unknown               Unknown  Unknown
libifcore.so.5     00002ABE6592EC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002ABE6592E14D  Unknown               Unknown  Unknown
libifcore.so.5     00002ABE65964340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002AFC7CD2FA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002AFC7CD2E4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002AFC7C0EE19E  Unknown               Unknown  Unknown
libifcore.so.5     00002AFC7C05DC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002AFC7C05D14D  Unknown               Unknown  Unknown
libifcore.so.5     00002AFC7C093340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B2934104A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B29341034B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B29334C319E  Unknown               Unknown  Unknown
libifcore.so.5     00002B2933432C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B293343214D  Unknown               Unknown  Unknown
libifcore.so.5     00002B2933468340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B3465EA8A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B3465EA74B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B346526719E  Unknown               Unknown  Unknown
libifcore.so.5     00002B34651D6C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B34651D614D  Unknown               Unknown  Unknown
libifcore.so.5     00002B346520C340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002ABD19699A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002ABD196984B6  Unknown               Unknown  Unknown
libifcore.so.5     00002ABD18A5819E  Unknown               Unknown  Unknown
libifcore.so.5     00002ABD189C7C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002ABD189C714D  Unknown               Unknown  Unknown
libifcore.so.5     00002ABD189FD340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B973D028A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B973D0274B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B973C3E719E  Unknown               Unknown  Unknown
libifcore.so.5     00002B973C356C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B973C35614D  Unknown               Unknown  Unknown
libifcore.so.5     00002B973C38C340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002AD72678EA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002AD72678D4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002AD725B4D19E  Unknown               Unknown  Unknown
libifcore.so.5     00002AD725ABCC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002AD725ABC14D  Unknown               Unknown  Unknown
libifcore.so.5     00002AD725AF2340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B611BBF6A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B611BBF54B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B611AFB519E  Unknown               Unknown  Unknown
libifcore.so.5     00002B611AF24C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B611AF2414D  Unknown               Unknown  Unknown
libifcore.so.5     00002B611AF5A340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B9257E2EA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B9257E2D4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B92571ED19E  Unknown               Unknown  Unknown
libifcore.so.5     00002B925715CC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B925715C14D  Unknown               Unknown  Unknown
libifcore.so.5     00002B9257192340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B719812FA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B719812E4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B71974EE19E  Unknown               Unknown  Unknown
libifcore.so.5     00002B719745DC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B719745D14D  Unknown               Unknown  Unknown
libifcore.so.5     00002B7197493340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B65F9E21A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B65F9E204B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B65F91E019E  Unknown               Unknown  Unknown
libifcore.so.5     00002B65F914FC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B65F914F14D  Unknown               Unknown  Unknown
libifcore.so.5     00002B65F9185340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B255F207A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B255F2064B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B255E5C619E  Unknown               Unknown  Unknown
libifcore.so.5     00002B255E535C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B255E53514D  Unknown               Unknown  Unknown
libifcore.so.5     00002B255E56B340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B5E25876A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B5E258754B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B5E24C3519E  Unknown               Unknown  Unknown
libifcore.so.5     00002B5E24BA4C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B5E24BA414D  Unknown               Unknown  Unknown
libifcore.so.5     00002B5E24BDA340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B11EDBEFA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B11EDBEE4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B11ECFAE19E  Unknown               Unknown  Unknown
libifcore.so.5     00002B11ECF1DC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B11ECF1D14D  Unknown               Unknown  Unknown
libifcore.so.5     00002B11ECF53340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
srun: error: curie1962: tasks 0-9,11-15: Exited with exit code 29
libintlc.so.5      00002B04A43FDA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B04A43FC4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B04A37BC19E  Unknown               Unknown  Unknown
libifcore.so.5     00002B04A372BC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B04A372B14D  Unknown               Unknown  Unknown
libifcore.so.5     00002B04A3761340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
srun: error: curie1964: tasks 16-17,19-22,24-28,30-31: Exited with exit code 29
srun: error: curie1962: task 10: Exited with exit code 29
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002BA96869EA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002BA96869D4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002BA967A5D19E  Unknown               Unknown  Unknown
libifcore.so.5     00002BA9679CCC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002BA9679CC14D  Unknown               Unknown  Unknown
libifcore.so.5     00002BA967A02340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
srun: error: curie1964: task 29: Exited with exit code 29
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002AB5C860FA1E  Unknown               Unknown  Unknown
libintlc.so.5      00002AB5C860E4B6  Unknown               Unknown  Unknown
libifcore.so.5     00002AB5C79CE19E  Unknown               Unknown  Unknown
libifcore.so.5     00002AB5C793DC4E  Unknown               Unknown  Unknown
libifcore.so.5     00002AB5C793D14D  Unknown               Unknown  Unknown
libifcore.so.5     00002AB5C7973340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 27, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/fort.27
Image              PC                Routine            Line        Source             
libintlc.so.5      00002B500E8D3A1E  Unknown               Unknown  Unknown
libintlc.so.5      00002B500E8D24B6  Unknown               Unknown  Unknown
libifcore.so.5     00002B500DC9219E  Unknown               Unknown  Unknown
libifcore.so.5     00002B500DC01C4E  Unknown               Unknown  Unknown
libifcore.so.5     00002B500DC0114D  Unknown               Unknown  Unknown
libifcore.so.5     00002B500DC37340  Unknown               Unknown  Unknown
xread_absorbing_i  0000000000408C7F  Unknown               Unknown  Unknown
srun: error: curie1964: tasks 18,23: Exited with exit code 29
++ date
+ echo Wed Jul 30 23:02:30 CEST 2014
+ echo ''
+ echo ' simulation'
++ date
+ echo Wed Jul 30 23:02:30 CEST 2014
+ run_simu
+ echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 3 !!!!!!!!!!!!!!!!'
+ pwd
+ ccc_mprun /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xspecfem3D
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000016_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002BA576A9ED1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000017_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B25221E6D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000018_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AB166A78D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000019_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B7C91471D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000021_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B9196968D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000022_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000000_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AE62913AD1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000023_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002ABE3B55BD1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AE9B5B48D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000001_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AE527B2FD1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000003_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B048E6EDD1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000004_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002ACB42846D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000006_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AAE8EE2AD1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000007_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AC2606C1D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000002_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B1F07457D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000009_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000026_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AE6CAEC2D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AE53555CD1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000012_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000024_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B6BD8460D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B0F8921ED1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000015_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000025_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B0BE9672D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B63665FCD1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000008_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000027_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AC30D188D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B41438F0D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000010_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000028_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002BA89C7ABD1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B682C318D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000011_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000029_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AAFFF4E1D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000030_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B70B1A95D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000031_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B58FEC33D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B4E89236D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000013_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AD9848E4D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000014_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B584F4FFD1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000005_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002B7294B80D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
srun: error: curie1964: tasks 16-19,21-31: Exited with exit code 47
srun: error: curie1962: tasks 0-15: Exited with exit code 47
forrtl: severe (47): write to READONLY file, unit 98, file /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/EXAMPLES/DSM_FOR_SPECFEM3D/example_simple_small/./OUTPUT_FILES/DATABASES_MPI/proc000020_res_minimum_period.vtk
Image              PC                Routine            Line        Source             
xspecfem3D         00000000006B061E  Unknown               Unknown  Unknown
xspecfem3D         00000000006AF0B6  Unknown               Unknown  Unknown
xspecfem3D         000000000066FCB2  Unknown               Unknown  Unknown
xspecfem3D         000000000061A58B  Unknown               Unknown  Unknown
xspecfem3D         0000000000619AF2  Unknown               Unknown  Unknown
xspecfem3D         000000000065BD2E  Unknown               Unknown  Unknown
xspecfem3D         0000000000602BD2  Unknown               Unknown  Unknown
xspecfem3D         00000000005C39E0  Unknown               Unknown  Unknown
xspecfem3D         0000000000568BF0  Unknown               Unknown  Unknown
xspecfem3D         0000000000596985  Unknown               Unknown  Unknown
xspecfem3D         0000000000561883  Unknown               Unknown  Unknown
xspecfem3D         000000000041C1BC  Unknown               Unknown  Unknown
libc.so.6          00002AFE48199D1D  Unknown               Unknown  Unknown
xspecfem3D         000000000041C0B9  Unknown               Unknown  Unknown
srun: error: curie1964: task 20: Exited with exit code 47
++ date
+ echo Wed Jul 30 23:02:35 CEST 2014
+ echo ''
+ echo ' MAKE movie'
++ date
+ echo Wed Jul 30 23:02:35 CEST 2014
+ create_movie velocity_Z_it /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 25 8500
+ cd bin
+ declare -i it itmax istep
+ PREFIX=velocity_Z_it
+ IN=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI
+ OUT=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie
+ istep=25
+ itmax=8500
+ mkdir /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie
mkdir: cannot create directory `/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie': File exists
+ it=25
+ '[' 25 -le 8500 ']'
+ '[' 25 -lt 1000000 ']'
+ FICHIER=velocity_Z_it25
+ '[' 25 -lt 100000 ']'
+ FICHIER=velocity_Z_it025
+ '[' 25 -lt 10000 ']'
+ FICHIER=velocity_Z_it0025
+ '[' 25 -lt 1000 ']'
+ FICHIER=velocity_Z_it00025
+ '[' 25 -lt 100 ']'
+ FICHIER=velocity_Z_it000025
+ echo velocity_Z_it000025
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000025 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=25+25
+ '[' 50 -le 8500 ']'
+ '[' 50 -lt 1000000 ']'
+ FICHIER=velocity_Z_it50
+ '[' 50 -lt 100000 ']'
+ FICHIER=velocity_Z_it050
+ '[' 50 -lt 10000 ']'
+ FICHIER=velocity_Z_it0050
+ '[' 50 -lt 1000 ']'
+ FICHIER=velocity_Z_it00050
+ '[' 50 -lt 100 ']'
+ FICHIER=velocity_Z_it000050
+ echo velocity_Z_it000050
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000050 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=50+25
+ '[' 75 -le 8500 ']'
+ '[' 75 -lt 1000000 ']'
+ FICHIER=velocity_Z_it75
+ '[' 75 -lt 100000 ']'
+ FICHIER=velocity_Z_it075
+ '[' 75 -lt 10000 ']'
+ FICHIER=velocity_Z_it0075
+ '[' 75 -lt 1000 ']'
+ FICHIER=velocity_Z_it00075
+ '[' 75 -lt 100 ']'
+ FICHIER=velocity_Z_it000075
+ echo velocity_Z_it000075
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000075 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=75+25
+ '[' 100 -le 8500 ']'
+ '[' 100 -lt 1000000 ']'
+ FICHIER=velocity_Z_it100
+ '[' 100 -lt 100000 ']'
+ FICHIER=velocity_Z_it0100
+ '[' 100 -lt 10000 ']'
+ FICHIER=velocity_Z_it00100
+ '[' 100 -lt 1000 ']'
+ FICHIER=velocity_Z_it000100
+ '[' 100 -lt 100 ']'
+ echo velocity_Z_it000100
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000100 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=100+25
+ '[' 125 -le 8500 ']'
+ '[' 125 -lt 1000000 ']'
+ FICHIER=velocity_Z_it125
+ '[' 125 -lt 100000 ']'
+ FICHIER=velocity_Z_it0125
+ '[' 125 -lt 10000 ']'
+ FICHIER=velocity_Z_it00125
+ '[' 125 -lt 1000 ']'
+ FICHIER=velocity_Z_it000125
+ '[' 125 -lt 100 ']'
+ echo velocity_Z_it000125
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000125 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=125+25
+ '[' 150 -le 8500 ']'
+ '[' 150 -lt 1000000 ']'
+ FICHIER=velocity_Z_it150
+ '[' 150 -lt 100000 ']'
+ FICHIER=velocity_Z_it0150
+ '[' 150 -lt 10000 ']'
+ FICHIER=velocity_Z_it00150
+ '[' 150 -lt 1000 ']'
+ FICHIER=velocity_Z_it000150
+ '[' 150 -lt 100 ']'
+ echo velocity_Z_it000150
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000150 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=150+25
+ '[' 175 -le 8500 ']'
+ '[' 175 -lt 1000000 ']'
+ FICHIER=velocity_Z_it175
+ '[' 175 -lt 100000 ']'
+ FICHIER=velocity_Z_it0175
+ '[' 175 -lt 10000 ']'
+ FICHIER=velocity_Z_it00175
+ '[' 175 -lt 1000 ']'
+ FICHIER=velocity_Z_it000175
+ '[' 175 -lt 100 ']'
+ echo velocity_Z_it000175
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000175 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=175+25
+ '[' 200 -le 8500 ']'
+ '[' 200 -lt 1000000 ']'
+ FICHIER=velocity_Z_it200
+ '[' 200 -lt 100000 ']'
+ FICHIER=velocity_Z_it0200
+ '[' 200 -lt 10000 ']'
+ FICHIER=velocity_Z_it00200
+ '[' 200 -lt 1000 ']'
+ FICHIER=velocity_Z_it000200
+ '[' 200 -lt 100 ']'
+ echo velocity_Z_it000200
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000200 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=200+25
+ '[' 225 -le 8500 ']'
+ '[' 225 -lt 1000000 ']'
+ FICHIER=velocity_Z_it225
+ '[' 225 -lt 100000 ']'
+ FICHIER=velocity_Z_it0225
+ '[' 225 -lt 10000 ']'
+ FICHIER=velocity_Z_it00225
+ '[' 225 -lt 1000 ']'
+ FICHIER=velocity_Z_it000225
+ '[' 225 -lt 100 ']'
+ echo velocity_Z_it000225
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000225 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=225+25
+ '[' 250 -le 8500 ']'
+ '[' 250 -lt 1000000 ']'
+ FICHIER=velocity_Z_it250
+ '[' 250 -lt 100000 ']'
+ FICHIER=velocity_Z_it0250
+ '[' 250 -lt 10000 ']'
+ FICHIER=velocity_Z_it00250
+ '[' 250 -lt 1000 ']'
+ FICHIER=velocity_Z_it000250
+ '[' 250 -lt 100 ']'
+ echo velocity_Z_it000250
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000250 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=250+25
+ '[' 275 -le 8500 ']'
+ '[' 275 -lt 1000000 ']'
+ FICHIER=velocity_Z_it275
+ '[' 275 -lt 100000 ']'
+ FICHIER=velocity_Z_it0275
+ '[' 275 -lt 10000 ']'
+ FICHIER=velocity_Z_it00275
+ '[' 275 -lt 1000 ']'
+ FICHIER=velocity_Z_it000275
+ '[' 275 -lt 100 ']'
+ echo velocity_Z_it000275
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000275 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=275+25
+ '[' 300 -le 8500 ']'
+ '[' 300 -lt 1000000 ']'
+ FICHIER=velocity_Z_it300
+ '[' 300 -lt 100000 ']'
+ FICHIER=velocity_Z_it0300
+ '[' 300 -lt 10000 ']'
+ FICHIER=velocity_Z_it00300
+ '[' 300 -lt 1000 ']'
+ FICHIER=velocity_Z_it000300
+ '[' 300 -lt 100 ']'
+ echo velocity_Z_it000300
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000300 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=300+25
+ '[' 325 -le 8500 ']'
+ '[' 325 -lt 1000000 ']'
+ FICHIER=velocity_Z_it325
+ '[' 325 -lt 100000 ']'
+ FICHIER=velocity_Z_it0325
+ '[' 325 -lt 10000 ']'
+ FICHIER=velocity_Z_it00325
+ '[' 325 -lt 1000 ']'
+ FICHIER=velocity_Z_it000325
+ '[' 325 -lt 100 ']'
+ echo velocity_Z_it000325
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000325 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=325+25
+ '[' 350 -le 8500 ']'
+ '[' 350 -lt 1000000 ']'
+ FICHIER=velocity_Z_it350
+ '[' 350 -lt 100000 ']'
+ FICHIER=velocity_Z_it0350
+ '[' 350 -lt 10000 ']'
+ FICHIER=velocity_Z_it00350
+ '[' 350 -lt 1000 ']'
+ FICHIER=velocity_Z_it000350
+ '[' 350 -lt 100 ']'
+ echo velocity_Z_it000350
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000350 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=350+25
+ '[' 375 -le 8500 ']'
+ '[' 375 -lt 1000000 ']'
+ FICHIER=velocity_Z_it375
+ '[' 375 -lt 100000 ']'
+ FICHIER=velocity_Z_it0375
+ '[' 375 -lt 10000 ']'
+ FICHIER=velocity_Z_it00375
+ '[' 375 -lt 1000 ']'
+ FICHIER=velocity_Z_it000375
+ '[' 375 -lt 100 ']'
+ echo velocity_Z_it000375
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000375 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=375+25
+ '[' 400 -le 8500 ']'
+ '[' 400 -lt 1000000 ']'
+ FICHIER=velocity_Z_it400
+ '[' 400 -lt 100000 ']'
+ FICHIER=velocity_Z_it0400
+ '[' 400 -lt 10000 ']'
+ FICHIER=velocity_Z_it00400
+ '[' 400 -lt 1000 ']'
+ FICHIER=velocity_Z_it000400
+ '[' 400 -lt 100 ']'
+ echo velocity_Z_it000400
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000400 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=400+25
+ '[' 425 -le 8500 ']'
+ '[' 425 -lt 1000000 ']'
+ FICHIER=velocity_Z_it425
+ '[' 425 -lt 100000 ']'
+ FICHIER=velocity_Z_it0425
+ '[' 425 -lt 10000 ']'
+ FICHIER=velocity_Z_it00425
+ '[' 425 -lt 1000 ']'
+ FICHIER=velocity_Z_it000425
+ '[' 425 -lt 100 ']'
+ echo velocity_Z_it000425
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000425 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=425+25
+ '[' 450 -le 8500 ']'
+ '[' 450 -lt 1000000 ']'
+ FICHIER=velocity_Z_it450
+ '[' 450 -lt 100000 ']'
+ FICHIER=velocity_Z_it0450
+ '[' 450 -lt 10000 ']'
+ FICHIER=velocity_Z_it00450
+ '[' 450 -lt 1000 ']'
+ FICHIER=velocity_Z_it000450
+ '[' 450 -lt 100 ']'
+ echo velocity_Z_it000450
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000450 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=450+25
+ '[' 475 -le 8500 ']'
+ '[' 475 -lt 1000000 ']'
+ FICHIER=velocity_Z_it475
+ '[' 475 -lt 100000 ']'
+ FICHIER=velocity_Z_it0475
+ '[' 475 -lt 10000 ']'
+ FICHIER=velocity_Z_it00475
+ '[' 475 -lt 1000 ']'
+ FICHIER=velocity_Z_it000475
+ '[' 475 -lt 100 ']'
+ echo velocity_Z_it000475
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000475 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=475+25
+ '[' 500 -le 8500 ']'
+ '[' 500 -lt 1000000 ']'
+ FICHIER=velocity_Z_it500
+ '[' 500 -lt 100000 ']'
+ FICHIER=velocity_Z_it0500
+ '[' 500 -lt 10000 ']'
+ FICHIER=velocity_Z_it00500
+ '[' 500 -lt 1000 ']'
+ FICHIER=velocity_Z_it000500
+ '[' 500 -lt 100 ']'
+ echo velocity_Z_it000500
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000500 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=500+25
+ '[' 525 -le 8500 ']'
+ '[' 525 -lt 1000000 ']'
+ FICHIER=velocity_Z_it525
+ '[' 525 -lt 100000 ']'
+ FICHIER=velocity_Z_it0525
+ '[' 525 -lt 10000 ']'
+ FICHIER=velocity_Z_it00525
+ '[' 525 -lt 1000 ']'
+ FICHIER=velocity_Z_it000525
+ '[' 525 -lt 100 ']'
+ echo velocity_Z_it000525
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000525 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=525+25
+ '[' 550 -le 8500 ']'
+ '[' 550 -lt 1000000 ']'
+ FICHIER=velocity_Z_it550
+ '[' 550 -lt 100000 ']'
+ FICHIER=velocity_Z_it0550
+ '[' 550 -lt 10000 ']'
+ FICHIER=velocity_Z_it00550
+ '[' 550 -lt 1000 ']'
+ FICHIER=velocity_Z_it000550
+ '[' 550 -lt 100 ']'
+ echo velocity_Z_it000550
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000550 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=550+25
+ '[' 575 -le 8500 ']'
+ '[' 575 -lt 1000000 ']'
+ FICHIER=velocity_Z_it575
+ '[' 575 -lt 100000 ']'
+ FICHIER=velocity_Z_it0575
+ '[' 575 -lt 10000 ']'
+ FICHIER=velocity_Z_it00575
+ '[' 575 -lt 1000 ']'
+ FICHIER=velocity_Z_it000575
+ '[' 575 -lt 100 ']'
+ echo velocity_Z_it000575
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000575 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=575+25
+ '[' 600 -le 8500 ']'
+ '[' 600 -lt 1000000 ']'
+ FICHIER=velocity_Z_it600
+ '[' 600 -lt 100000 ']'
+ FICHIER=velocity_Z_it0600
+ '[' 600 -lt 10000 ']'
+ FICHIER=velocity_Z_it00600
+ '[' 600 -lt 1000 ']'
+ FICHIER=velocity_Z_it000600
+ '[' 600 -lt 100 ']'
+ echo velocity_Z_it000600
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000600 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=600+25
+ '[' 625 -le 8500 ']'
+ '[' 625 -lt 1000000 ']'
+ FICHIER=velocity_Z_it625
+ '[' 625 -lt 100000 ']'
+ FICHIER=velocity_Z_it0625
+ '[' 625 -lt 10000 ']'
+ FICHIER=velocity_Z_it00625
+ '[' 625 -lt 1000 ']'
+ FICHIER=velocity_Z_it000625
+ '[' 625 -lt 100 ']'
+ echo velocity_Z_it000625
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000625 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=625+25
+ '[' 650 -le 8500 ']'
+ '[' 650 -lt 1000000 ']'
+ FICHIER=velocity_Z_it650
+ '[' 650 -lt 100000 ']'
+ FICHIER=velocity_Z_it0650
+ '[' 650 -lt 10000 ']'
+ FICHIER=velocity_Z_it00650
+ '[' 650 -lt 1000 ']'
+ FICHIER=velocity_Z_it000650
+ '[' 650 -lt 100 ']'
+ echo velocity_Z_it000650
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000650 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=650+25
+ '[' 675 -le 8500 ']'
+ '[' 675 -lt 1000000 ']'
+ FICHIER=velocity_Z_it675
+ '[' 675 -lt 100000 ']'
+ FICHIER=velocity_Z_it0675
+ '[' 675 -lt 10000 ']'
+ FICHIER=velocity_Z_it00675
+ '[' 675 -lt 1000 ']'
+ FICHIER=velocity_Z_it000675
+ '[' 675 -lt 100 ']'
+ echo velocity_Z_it000675
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000675 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=675+25
+ '[' 700 -le 8500 ']'
+ '[' 700 -lt 1000000 ']'
+ FICHIER=velocity_Z_it700
+ '[' 700 -lt 100000 ']'
+ FICHIER=velocity_Z_it0700
+ '[' 700 -lt 10000 ']'
+ FICHIER=velocity_Z_it00700
+ '[' 700 -lt 1000 ']'
+ FICHIER=velocity_Z_it000700
+ '[' 700 -lt 100 ']'
+ echo velocity_Z_it000700
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000700 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=700+25
+ '[' 725 -le 8500 ']'
+ '[' 725 -lt 1000000 ']'
+ FICHIER=velocity_Z_it725
+ '[' 725 -lt 100000 ']'
+ FICHIER=velocity_Z_it0725
+ '[' 725 -lt 10000 ']'
+ FICHIER=velocity_Z_it00725
+ '[' 725 -lt 1000 ']'
+ FICHIER=velocity_Z_it000725
+ '[' 725 -lt 100 ']'
+ echo velocity_Z_it000725
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000725 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=725+25
+ '[' 750 -le 8500 ']'
+ '[' 750 -lt 1000000 ']'
+ FICHIER=velocity_Z_it750
+ '[' 750 -lt 100000 ']'
+ FICHIER=velocity_Z_it0750
+ '[' 750 -lt 10000 ']'
+ FICHIER=velocity_Z_it00750
+ '[' 750 -lt 1000 ']'
+ FICHIER=velocity_Z_it000750
+ '[' 750 -lt 100 ']'
+ echo velocity_Z_it000750
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000750 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=750+25
+ '[' 775 -le 8500 ']'
+ '[' 775 -lt 1000000 ']'
+ FICHIER=velocity_Z_it775
+ '[' 775 -lt 100000 ']'
+ FICHIER=velocity_Z_it0775
+ '[' 775 -lt 10000 ']'
+ FICHIER=velocity_Z_it00775
+ '[' 775 -lt 1000 ']'
+ FICHIER=velocity_Z_it000775
+ '[' 775 -lt 100 ']'
+ echo velocity_Z_it000775
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000775 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=775+25
+ '[' 800 -le 8500 ']'
+ '[' 800 -lt 1000000 ']'
+ FICHIER=velocity_Z_it800
+ '[' 800 -lt 100000 ']'
+ FICHIER=velocity_Z_it0800
+ '[' 800 -lt 10000 ']'
+ FICHIER=velocity_Z_it00800
+ '[' 800 -lt 1000 ']'
+ FICHIER=velocity_Z_it000800
+ '[' 800 -lt 100 ']'
+ echo velocity_Z_it000800
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000800 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=800+25
+ '[' 825 -le 8500 ']'
+ '[' 825 -lt 1000000 ']'
+ FICHIER=velocity_Z_it825
+ '[' 825 -lt 100000 ']'
+ FICHIER=velocity_Z_it0825
+ '[' 825 -lt 10000 ']'
+ FICHIER=velocity_Z_it00825
+ '[' 825 -lt 1000 ']'
+ FICHIER=velocity_Z_it000825
+ '[' 825 -lt 100 ']'
+ echo velocity_Z_it000825
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000825 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=825+25
+ '[' 850 -le 8500 ']'
+ '[' 850 -lt 1000000 ']'
+ FICHIER=velocity_Z_it850
+ '[' 850 -lt 100000 ']'
+ FICHIER=velocity_Z_it0850
+ '[' 850 -lt 10000 ']'
+ FICHIER=velocity_Z_it00850
+ '[' 850 -lt 1000 ']'
+ FICHIER=velocity_Z_it000850
+ '[' 850 -lt 100 ']'
+ echo velocity_Z_it000850
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000850 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=850+25
+ '[' 875 -le 8500 ']'
+ '[' 875 -lt 1000000 ']'
+ FICHIER=velocity_Z_it875
+ '[' 875 -lt 100000 ']'
+ FICHIER=velocity_Z_it0875
+ '[' 875 -lt 10000 ']'
+ FICHIER=velocity_Z_it00875
+ '[' 875 -lt 1000 ']'
+ FICHIER=velocity_Z_it000875
+ '[' 875 -lt 100 ']'
+ echo velocity_Z_it000875
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000875 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=875+25
+ '[' 900 -le 8500 ']'
+ '[' 900 -lt 1000000 ']'
+ FICHIER=velocity_Z_it900
+ '[' 900 -lt 100000 ']'
+ FICHIER=velocity_Z_it0900
+ '[' 900 -lt 10000 ']'
+ FICHIER=velocity_Z_it00900
+ '[' 900 -lt 1000 ']'
+ FICHIER=velocity_Z_it000900
+ '[' 900 -lt 100 ']'
+ echo velocity_Z_it000900
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000900 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=900+25
+ '[' 925 -le 8500 ']'
+ '[' 925 -lt 1000000 ']'
+ FICHIER=velocity_Z_it925
+ '[' 925 -lt 100000 ']'
+ FICHIER=velocity_Z_it0925
+ '[' 925 -lt 10000 ']'
+ FICHIER=velocity_Z_it00925
+ '[' 925 -lt 1000 ']'
+ FICHIER=velocity_Z_it000925
+ '[' 925 -lt 100 ']'
+ echo velocity_Z_it000925
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000925 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=925+25
+ '[' 950 -le 8500 ']'
+ '[' 950 -lt 1000000 ']'
+ FICHIER=velocity_Z_it950
+ '[' 950 -lt 100000 ']'
+ FICHIER=velocity_Z_it0950
+ '[' 950 -lt 10000 ']'
+ FICHIER=velocity_Z_it00950
+ '[' 950 -lt 1000 ']'
+ FICHIER=velocity_Z_it000950
+ '[' 950 -lt 100 ']'
+ echo velocity_Z_it000950
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000950 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=950+25
+ '[' 975 -le 8500 ']'
+ '[' 975 -lt 1000000 ']'
+ FICHIER=velocity_Z_it975
+ '[' 975 -lt 100000 ']'
+ FICHIER=velocity_Z_it0975
+ '[' 975 -lt 10000 ']'
+ FICHIER=velocity_Z_it00975
+ '[' 975 -lt 1000 ']'
+ FICHIER=velocity_Z_it000975
+ '[' 975 -lt 100 ']'
+ echo velocity_Z_it000975
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it000975 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=975+25
+ '[' 1000 -le 8500 ']'
+ '[' 1000 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1000
+ '[' 1000 -lt 100000 ']'
+ FICHIER=velocity_Z_it01000
+ '[' 1000 -lt 10000 ']'
+ FICHIER=velocity_Z_it001000
+ '[' 1000 -lt 1000 ']'
+ '[' 1000 -lt 100 ']'
+ echo velocity_Z_it001000
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001000 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1000+25
+ '[' 1025 -le 8500 ']'
+ '[' 1025 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1025
+ '[' 1025 -lt 100000 ']'
+ FICHIER=velocity_Z_it01025
+ '[' 1025 -lt 10000 ']'
+ FICHIER=velocity_Z_it001025
+ '[' 1025 -lt 1000 ']'
+ '[' 1025 -lt 100 ']'
+ echo velocity_Z_it001025
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001025 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1025+25
+ '[' 1050 -le 8500 ']'
+ '[' 1050 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1050
+ '[' 1050 -lt 100000 ']'
+ FICHIER=velocity_Z_it01050
+ '[' 1050 -lt 10000 ']'
+ FICHIER=velocity_Z_it001050
+ '[' 1050 -lt 1000 ']'
+ '[' 1050 -lt 100 ']'
+ echo velocity_Z_it001050
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001050 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1050+25
+ '[' 1075 -le 8500 ']'
+ '[' 1075 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1075
+ '[' 1075 -lt 100000 ']'
+ FICHIER=velocity_Z_it01075
+ '[' 1075 -lt 10000 ']'
+ FICHIER=velocity_Z_it001075
+ '[' 1075 -lt 1000 ']'
+ '[' 1075 -lt 100 ']'
+ echo velocity_Z_it001075
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001075 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1075+25
+ '[' 1100 -le 8500 ']'
+ '[' 1100 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1100
+ '[' 1100 -lt 100000 ']'
+ FICHIER=velocity_Z_it01100
+ '[' 1100 -lt 10000 ']'
+ FICHIER=velocity_Z_it001100
+ '[' 1100 -lt 1000 ']'
+ '[' 1100 -lt 100 ']'
+ echo velocity_Z_it001100
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001100 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1100+25
+ '[' 1125 -le 8500 ']'
+ '[' 1125 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1125
+ '[' 1125 -lt 100000 ']'
+ FICHIER=velocity_Z_it01125
+ '[' 1125 -lt 10000 ']'
+ FICHIER=velocity_Z_it001125
+ '[' 1125 -lt 1000 ']'
+ '[' 1125 -lt 100 ']'
+ echo velocity_Z_it001125
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001125 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1125+25
+ '[' 1150 -le 8500 ']'
+ '[' 1150 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1150
+ '[' 1150 -lt 100000 ']'
+ FICHIER=velocity_Z_it01150
+ '[' 1150 -lt 10000 ']'
+ FICHIER=velocity_Z_it001150
+ '[' 1150 -lt 1000 ']'
+ '[' 1150 -lt 100 ']'
+ echo velocity_Z_it001150
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001150 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1150+25
+ '[' 1175 -le 8500 ']'
+ '[' 1175 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1175
+ '[' 1175 -lt 100000 ']'
+ FICHIER=velocity_Z_it01175
+ '[' 1175 -lt 10000 ']'
+ FICHIER=velocity_Z_it001175
+ '[' 1175 -lt 1000 ']'
+ '[' 1175 -lt 100 ']'
+ echo velocity_Z_it001175
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001175 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1175+25
+ '[' 1200 -le 8500 ']'
+ '[' 1200 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1200
+ '[' 1200 -lt 100000 ']'
+ FICHIER=velocity_Z_it01200
+ '[' 1200 -lt 10000 ']'
+ FICHIER=velocity_Z_it001200
+ '[' 1200 -lt 1000 ']'
+ '[' 1200 -lt 100 ']'
+ echo velocity_Z_it001200
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001200 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1200+25
+ '[' 1225 -le 8500 ']'
+ '[' 1225 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1225
+ '[' 1225 -lt 100000 ']'
+ FICHIER=velocity_Z_it01225
+ '[' 1225 -lt 10000 ']'
+ FICHIER=velocity_Z_it001225
+ '[' 1225 -lt 1000 ']'
+ '[' 1225 -lt 100 ']'
+ echo velocity_Z_it001225
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001225 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1225+25
+ '[' 1250 -le 8500 ']'
+ '[' 1250 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1250
+ '[' 1250 -lt 100000 ']'
+ FICHIER=velocity_Z_it01250
+ '[' 1250 -lt 10000 ']'
+ FICHIER=velocity_Z_it001250
+ '[' 1250 -lt 1000 ']'
+ '[' 1250 -lt 100 ']'
+ echo velocity_Z_it001250
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001250 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1250+25
+ '[' 1275 -le 8500 ']'
+ '[' 1275 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1275
+ '[' 1275 -lt 100000 ']'
+ FICHIER=velocity_Z_it01275
+ '[' 1275 -lt 10000 ']'
+ FICHIER=velocity_Z_it001275
+ '[' 1275 -lt 1000 ']'
+ '[' 1275 -lt 100 ']'
+ echo velocity_Z_it001275
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001275 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1275+25
+ '[' 1300 -le 8500 ']'
+ '[' 1300 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1300
+ '[' 1300 -lt 100000 ']'
+ FICHIER=velocity_Z_it01300
+ '[' 1300 -lt 10000 ']'
+ FICHIER=velocity_Z_it001300
+ '[' 1300 -lt 1000 ']'
+ '[' 1300 -lt 100 ']'
+ echo velocity_Z_it001300
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001300 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1300+25
+ '[' 1325 -le 8500 ']'
+ '[' 1325 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1325
+ '[' 1325 -lt 100000 ']'
+ FICHIER=velocity_Z_it01325
+ '[' 1325 -lt 10000 ']'
+ FICHIER=velocity_Z_it001325
+ '[' 1325 -lt 1000 ']'
+ '[' 1325 -lt 100 ']'
+ echo velocity_Z_it001325
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001325 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1325+25
+ '[' 1350 -le 8500 ']'
+ '[' 1350 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1350
+ '[' 1350 -lt 100000 ']'
+ FICHIER=velocity_Z_it01350
+ '[' 1350 -lt 10000 ']'
+ FICHIER=velocity_Z_it001350
+ '[' 1350 -lt 1000 ']'
+ '[' 1350 -lt 100 ']'
+ echo velocity_Z_it001350
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001350 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1350+25
+ '[' 1375 -le 8500 ']'
+ '[' 1375 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1375
+ '[' 1375 -lt 100000 ']'
+ FICHIER=velocity_Z_it01375
+ '[' 1375 -lt 10000 ']'
+ FICHIER=velocity_Z_it001375
+ '[' 1375 -lt 1000 ']'
+ '[' 1375 -lt 100 ']'
+ echo velocity_Z_it001375
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001375 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1375+25
+ '[' 1400 -le 8500 ']'
+ '[' 1400 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1400
+ '[' 1400 -lt 100000 ']'
+ FICHIER=velocity_Z_it01400
+ '[' 1400 -lt 10000 ']'
+ FICHIER=velocity_Z_it001400
+ '[' 1400 -lt 1000 ']'
+ '[' 1400 -lt 100 ']'
+ echo velocity_Z_it001400
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001400 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1400+25
+ '[' 1425 -le 8500 ']'
+ '[' 1425 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1425
+ '[' 1425 -lt 100000 ']'
+ FICHIER=velocity_Z_it01425
+ '[' 1425 -lt 10000 ']'
+ FICHIER=velocity_Z_it001425
+ '[' 1425 -lt 1000 ']'
+ '[' 1425 -lt 100 ']'
+ echo velocity_Z_it001425
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001425 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1425+25
+ '[' 1450 -le 8500 ']'
+ '[' 1450 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1450
+ '[' 1450 -lt 100000 ']'
+ FICHIER=velocity_Z_it01450
+ '[' 1450 -lt 10000 ']'
+ FICHIER=velocity_Z_it001450
+ '[' 1450 -lt 1000 ']'
+ '[' 1450 -lt 100 ']'
+ echo velocity_Z_it001450
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001450 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1450+25
+ '[' 1475 -le 8500 ']'
+ '[' 1475 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1475
+ '[' 1475 -lt 100000 ']'
+ FICHIER=velocity_Z_it01475
+ '[' 1475 -lt 10000 ']'
+ FICHIER=velocity_Z_it001475
+ '[' 1475 -lt 1000 ']'
+ '[' 1475 -lt 100 ']'
+ echo velocity_Z_it001475
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001475 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1475+25
+ '[' 1500 -le 8500 ']'
+ '[' 1500 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1500
+ '[' 1500 -lt 100000 ']'
+ FICHIER=velocity_Z_it01500
+ '[' 1500 -lt 10000 ']'
+ FICHIER=velocity_Z_it001500
+ '[' 1500 -lt 1000 ']'
+ '[' 1500 -lt 100 ']'
+ echo velocity_Z_it001500
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001500 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1500+25
+ '[' 1525 -le 8500 ']'
+ '[' 1525 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1525
+ '[' 1525 -lt 100000 ']'
+ FICHIER=velocity_Z_it01525
+ '[' 1525 -lt 10000 ']'
+ FICHIER=velocity_Z_it001525
+ '[' 1525 -lt 1000 ']'
+ '[' 1525 -lt 100 ']'
+ echo velocity_Z_it001525
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001525 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1525+25
+ '[' 1550 -le 8500 ']'
+ '[' 1550 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1550
+ '[' 1550 -lt 100000 ']'
+ FICHIER=velocity_Z_it01550
+ '[' 1550 -lt 10000 ']'
+ FICHIER=velocity_Z_it001550
+ '[' 1550 -lt 1000 ']'
+ '[' 1550 -lt 100 ']'
+ echo velocity_Z_it001550
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001550 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1550+25
+ '[' 1575 -le 8500 ']'
+ '[' 1575 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1575
+ '[' 1575 -lt 100000 ']'
+ FICHIER=velocity_Z_it01575
+ '[' 1575 -lt 10000 ']'
+ FICHIER=velocity_Z_it001575
+ '[' 1575 -lt 1000 ']'
+ '[' 1575 -lt 100 ']'
+ echo velocity_Z_it001575
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001575 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1575+25
+ '[' 1600 -le 8500 ']'
+ '[' 1600 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1600
+ '[' 1600 -lt 100000 ']'
+ FICHIER=velocity_Z_it01600
+ '[' 1600 -lt 10000 ']'
+ FICHIER=velocity_Z_it001600
+ '[' 1600 -lt 1000 ']'
+ '[' 1600 -lt 100 ']'
+ echo velocity_Z_it001600
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001600 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1600+25
+ '[' 1625 -le 8500 ']'
+ '[' 1625 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1625
+ '[' 1625 -lt 100000 ']'
+ FICHIER=velocity_Z_it01625
+ '[' 1625 -lt 10000 ']'
+ FICHIER=velocity_Z_it001625
+ '[' 1625 -lt 1000 ']'
+ '[' 1625 -lt 100 ']'
+ echo velocity_Z_it001625
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001625 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1625+25
+ '[' 1650 -le 8500 ']'
+ '[' 1650 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1650
+ '[' 1650 -lt 100000 ']'
+ FICHIER=velocity_Z_it01650
+ '[' 1650 -lt 10000 ']'
+ FICHIER=velocity_Z_it001650
+ '[' 1650 -lt 1000 ']'
+ '[' 1650 -lt 100 ']'
+ echo velocity_Z_it001650
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001650 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1650+25
+ '[' 1675 -le 8500 ']'
+ '[' 1675 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1675
+ '[' 1675 -lt 100000 ']'
+ FICHIER=velocity_Z_it01675
+ '[' 1675 -lt 10000 ']'
+ FICHIER=velocity_Z_it001675
+ '[' 1675 -lt 1000 ']'
+ '[' 1675 -lt 100 ']'
+ echo velocity_Z_it001675
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001675 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1675+25
+ '[' 1700 -le 8500 ']'
+ '[' 1700 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1700
+ '[' 1700 -lt 100000 ']'
+ FICHIER=velocity_Z_it01700
+ '[' 1700 -lt 10000 ']'
+ FICHIER=velocity_Z_it001700
+ '[' 1700 -lt 1000 ']'
+ '[' 1700 -lt 100 ']'
+ echo velocity_Z_it001700
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001700 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1700+25
+ '[' 1725 -le 8500 ']'
+ '[' 1725 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1725
+ '[' 1725 -lt 100000 ']'
+ FICHIER=velocity_Z_it01725
+ '[' 1725 -lt 10000 ']'
+ FICHIER=velocity_Z_it001725
+ '[' 1725 -lt 1000 ']'
+ '[' 1725 -lt 100 ']'
+ echo velocity_Z_it001725
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001725 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1725+25
+ '[' 1750 -le 8500 ']'
+ '[' 1750 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1750
+ '[' 1750 -lt 100000 ']'
+ FICHIER=velocity_Z_it01750
+ '[' 1750 -lt 10000 ']'
+ FICHIER=velocity_Z_it001750
+ '[' 1750 -lt 1000 ']'
+ '[' 1750 -lt 100 ']'
+ echo velocity_Z_it001750
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001750 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1750+25
+ '[' 1775 -le 8500 ']'
+ '[' 1775 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1775
+ '[' 1775 -lt 100000 ']'
+ FICHIER=velocity_Z_it01775
+ '[' 1775 -lt 10000 ']'
+ FICHIER=velocity_Z_it001775
+ '[' 1775 -lt 1000 ']'
+ '[' 1775 -lt 100 ']'
+ echo velocity_Z_it001775
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001775 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1775+25
+ '[' 1800 -le 8500 ']'
+ '[' 1800 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1800
+ '[' 1800 -lt 100000 ']'
+ FICHIER=velocity_Z_it01800
+ '[' 1800 -lt 10000 ']'
+ FICHIER=velocity_Z_it001800
+ '[' 1800 -lt 1000 ']'
+ '[' 1800 -lt 100 ']'
+ echo velocity_Z_it001800
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001800 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1800+25
+ '[' 1825 -le 8500 ']'
+ '[' 1825 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1825
+ '[' 1825 -lt 100000 ']'
+ FICHIER=velocity_Z_it01825
+ '[' 1825 -lt 10000 ']'
+ FICHIER=velocity_Z_it001825
+ '[' 1825 -lt 1000 ']'
+ '[' 1825 -lt 100 ']'
+ echo velocity_Z_it001825
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001825 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1825+25
+ '[' 1850 -le 8500 ']'
+ '[' 1850 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1850
+ '[' 1850 -lt 100000 ']'
+ FICHIER=velocity_Z_it01850
+ '[' 1850 -lt 10000 ']'
+ FICHIER=velocity_Z_it001850
+ '[' 1850 -lt 1000 ']'
+ '[' 1850 -lt 100 ']'
+ echo velocity_Z_it001850
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001850 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1850+25
+ '[' 1875 -le 8500 ']'
+ '[' 1875 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1875
+ '[' 1875 -lt 100000 ']'
+ FICHIER=velocity_Z_it01875
+ '[' 1875 -lt 10000 ']'
+ FICHIER=velocity_Z_it001875
+ '[' 1875 -lt 1000 ']'
+ '[' 1875 -lt 100 ']'
+ echo velocity_Z_it001875
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001875 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1875+25
+ '[' 1900 -le 8500 ']'
+ '[' 1900 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1900
+ '[' 1900 -lt 100000 ']'
+ FICHIER=velocity_Z_it01900
+ '[' 1900 -lt 10000 ']'
+ FICHIER=velocity_Z_it001900
+ '[' 1900 -lt 1000 ']'
+ '[' 1900 -lt 100 ']'
+ echo velocity_Z_it001900
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001900 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1900+25
+ '[' 1925 -le 8500 ']'
+ '[' 1925 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1925
+ '[' 1925 -lt 100000 ']'
+ FICHIER=velocity_Z_it01925
+ '[' 1925 -lt 10000 ']'
+ FICHIER=velocity_Z_it001925
+ '[' 1925 -lt 1000 ']'
+ '[' 1925 -lt 100 ']'
+ echo velocity_Z_it001925
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001925 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1925+25
+ '[' 1950 -le 8500 ']'
+ '[' 1950 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1950
+ '[' 1950 -lt 100000 ']'
+ FICHIER=velocity_Z_it01950
+ '[' 1950 -lt 10000 ']'
+ FICHIER=velocity_Z_it001950
+ '[' 1950 -lt 1000 ']'
+ '[' 1950 -lt 100 ']'
+ echo velocity_Z_it001950
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001950 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1950+25
+ '[' 1975 -le 8500 ']'
+ '[' 1975 -lt 1000000 ']'
+ FICHIER=velocity_Z_it1975
+ '[' 1975 -lt 100000 ']'
+ FICHIER=velocity_Z_it01975
+ '[' 1975 -lt 10000 ']'
+ FICHIER=velocity_Z_it001975
+ '[' 1975 -lt 1000 ']'
+ '[' 1975 -lt 100 ']'
+ echo velocity_Z_it001975
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it001975 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=1975+25
+ '[' 2000 -le 8500 ']'
+ '[' 2000 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2000
+ '[' 2000 -lt 100000 ']'
+ FICHIER=velocity_Z_it02000
+ '[' 2000 -lt 10000 ']'
+ FICHIER=velocity_Z_it002000
+ '[' 2000 -lt 1000 ']'
+ '[' 2000 -lt 100 ']'
+ echo velocity_Z_it002000
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002000 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2000+25
+ '[' 2025 -le 8500 ']'
+ '[' 2025 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2025
+ '[' 2025 -lt 100000 ']'
+ FICHIER=velocity_Z_it02025
+ '[' 2025 -lt 10000 ']'
+ FICHIER=velocity_Z_it002025
+ '[' 2025 -lt 1000 ']'
+ '[' 2025 -lt 100 ']'
+ echo velocity_Z_it002025
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002025 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2025+25
+ '[' 2050 -le 8500 ']'
+ '[' 2050 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2050
+ '[' 2050 -lt 100000 ']'
+ FICHIER=velocity_Z_it02050
+ '[' 2050 -lt 10000 ']'
+ FICHIER=velocity_Z_it002050
+ '[' 2050 -lt 1000 ']'
+ '[' 2050 -lt 100 ']'
+ echo velocity_Z_it002050
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002050 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2050+25
+ '[' 2075 -le 8500 ']'
+ '[' 2075 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2075
+ '[' 2075 -lt 100000 ']'
+ FICHIER=velocity_Z_it02075
+ '[' 2075 -lt 10000 ']'
+ FICHIER=velocity_Z_it002075
+ '[' 2075 -lt 1000 ']'
+ '[' 2075 -lt 100 ']'
+ echo velocity_Z_it002075
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002075 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2075+25
+ '[' 2100 -le 8500 ']'
+ '[' 2100 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2100
+ '[' 2100 -lt 100000 ']'
+ FICHIER=velocity_Z_it02100
+ '[' 2100 -lt 10000 ']'
+ FICHIER=velocity_Z_it002100
+ '[' 2100 -lt 1000 ']'
+ '[' 2100 -lt 100 ']'
+ echo velocity_Z_it002100
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002100 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2100+25
+ '[' 2125 -le 8500 ']'
+ '[' 2125 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2125
+ '[' 2125 -lt 100000 ']'
+ FICHIER=velocity_Z_it02125
+ '[' 2125 -lt 10000 ']'
+ FICHIER=velocity_Z_it002125
+ '[' 2125 -lt 1000 ']'
+ '[' 2125 -lt 100 ']'
+ echo velocity_Z_it002125
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002125 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2125+25
+ '[' 2150 -le 8500 ']'
+ '[' 2150 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2150
+ '[' 2150 -lt 100000 ']'
+ FICHIER=velocity_Z_it02150
+ '[' 2150 -lt 10000 ']'
+ FICHIER=velocity_Z_it002150
+ '[' 2150 -lt 1000 ']'
+ '[' 2150 -lt 100 ']'
+ echo velocity_Z_it002150
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002150 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2150+25
+ '[' 2175 -le 8500 ']'
+ '[' 2175 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2175
+ '[' 2175 -lt 100000 ']'
+ FICHIER=velocity_Z_it02175
+ '[' 2175 -lt 10000 ']'
+ FICHIER=velocity_Z_it002175
+ '[' 2175 -lt 1000 ']'
+ '[' 2175 -lt 100 ']'
+ echo velocity_Z_it002175
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002175 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2175+25
+ '[' 2200 -le 8500 ']'
+ '[' 2200 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2200
+ '[' 2200 -lt 100000 ']'
+ FICHIER=velocity_Z_it02200
+ '[' 2200 -lt 10000 ']'
+ FICHIER=velocity_Z_it002200
+ '[' 2200 -lt 1000 ']'
+ '[' 2200 -lt 100 ']'
+ echo velocity_Z_it002200
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002200 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2200+25
+ '[' 2225 -le 8500 ']'
+ '[' 2225 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2225
+ '[' 2225 -lt 100000 ']'
+ FICHIER=velocity_Z_it02225
+ '[' 2225 -lt 10000 ']'
+ FICHIER=velocity_Z_it002225
+ '[' 2225 -lt 1000 ']'
+ '[' 2225 -lt 100 ']'
+ echo velocity_Z_it002225
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002225 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2225+25
+ '[' 2250 -le 8500 ']'
+ '[' 2250 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2250
+ '[' 2250 -lt 100000 ']'
+ FICHIER=velocity_Z_it02250
+ '[' 2250 -lt 10000 ']'
+ FICHIER=velocity_Z_it002250
+ '[' 2250 -lt 1000 ']'
+ '[' 2250 -lt 100 ']'
+ echo velocity_Z_it002250
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002250 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2250+25
+ '[' 2275 -le 8500 ']'
+ '[' 2275 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2275
+ '[' 2275 -lt 100000 ']'
+ FICHIER=velocity_Z_it02275
+ '[' 2275 -lt 10000 ']'
+ FICHIER=velocity_Z_it002275
+ '[' 2275 -lt 1000 ']'
+ '[' 2275 -lt 100 ']'
+ echo velocity_Z_it002275
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002275 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2275+25
+ '[' 2300 -le 8500 ']'
+ '[' 2300 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2300
+ '[' 2300 -lt 100000 ']'
+ FICHIER=velocity_Z_it02300
+ '[' 2300 -lt 10000 ']'
+ FICHIER=velocity_Z_it002300
+ '[' 2300 -lt 1000 ']'
+ '[' 2300 -lt 100 ']'
+ echo velocity_Z_it002300
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002300 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2300+25
+ '[' 2325 -le 8500 ']'
+ '[' 2325 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2325
+ '[' 2325 -lt 100000 ']'
+ FICHIER=velocity_Z_it02325
+ '[' 2325 -lt 10000 ']'
+ FICHIER=velocity_Z_it002325
+ '[' 2325 -lt 1000 ']'
+ '[' 2325 -lt 100 ']'
+ echo velocity_Z_it002325
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002325 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2325+25
+ '[' 2350 -le 8500 ']'
+ '[' 2350 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2350
+ '[' 2350 -lt 100000 ']'
+ FICHIER=velocity_Z_it02350
+ '[' 2350 -lt 10000 ']'
+ FICHIER=velocity_Z_it002350
+ '[' 2350 -lt 1000 ']'
+ '[' 2350 -lt 100 ']'
+ echo velocity_Z_it002350
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002350 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2350+25
+ '[' 2375 -le 8500 ']'
+ '[' 2375 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2375
+ '[' 2375 -lt 100000 ']'
+ FICHIER=velocity_Z_it02375
+ '[' 2375 -lt 10000 ']'
+ FICHIER=velocity_Z_it002375
+ '[' 2375 -lt 1000 ']'
+ '[' 2375 -lt 100 ']'
+ echo velocity_Z_it002375
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002375 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2375+25
+ '[' 2400 -le 8500 ']'
+ '[' 2400 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2400
+ '[' 2400 -lt 100000 ']'
+ FICHIER=velocity_Z_it02400
+ '[' 2400 -lt 10000 ']'
+ FICHIER=velocity_Z_it002400
+ '[' 2400 -lt 1000 ']'
+ '[' 2400 -lt 100 ']'
+ echo velocity_Z_it002400
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002400 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2400+25
+ '[' 2425 -le 8500 ']'
+ '[' 2425 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2425
+ '[' 2425 -lt 100000 ']'
+ FICHIER=velocity_Z_it02425
+ '[' 2425 -lt 10000 ']'
+ FICHIER=velocity_Z_it002425
+ '[' 2425 -lt 1000 ']'
+ '[' 2425 -lt 100 ']'
+ echo velocity_Z_it002425
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002425 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2425+25
+ '[' 2450 -le 8500 ']'
+ '[' 2450 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2450
+ '[' 2450 -lt 100000 ']'
+ FICHIER=velocity_Z_it02450
+ '[' 2450 -lt 10000 ']'
+ FICHIER=velocity_Z_it002450
+ '[' 2450 -lt 1000 ']'
+ '[' 2450 -lt 100 ']'
+ echo velocity_Z_it002450
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002450 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2450+25
+ '[' 2475 -le 8500 ']'
+ '[' 2475 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2475
+ '[' 2475 -lt 100000 ']'
+ FICHIER=velocity_Z_it02475
+ '[' 2475 -lt 10000 ']'
+ FICHIER=velocity_Z_it002475
+ '[' 2475 -lt 1000 ']'
+ '[' 2475 -lt 100 ']'
+ echo velocity_Z_it002475
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002475 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2475+25
+ '[' 2500 -le 8500 ']'
+ '[' 2500 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2500
+ '[' 2500 -lt 100000 ']'
+ FICHIER=velocity_Z_it02500
+ '[' 2500 -lt 10000 ']'
+ FICHIER=velocity_Z_it002500
+ '[' 2500 -lt 1000 ']'
+ '[' 2500 -lt 100 ']'
+ echo velocity_Z_it002500
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002500 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2500+25
+ '[' 2525 -le 8500 ']'
+ '[' 2525 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2525
+ '[' 2525 -lt 100000 ']'
+ FICHIER=velocity_Z_it02525
+ '[' 2525 -lt 10000 ']'
+ FICHIER=velocity_Z_it002525
+ '[' 2525 -lt 1000 ']'
+ '[' 2525 -lt 100 ']'
+ echo velocity_Z_it002525
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002525 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2525+25
+ '[' 2550 -le 8500 ']'
+ '[' 2550 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2550
+ '[' 2550 -lt 100000 ']'
+ FICHIER=velocity_Z_it02550
+ '[' 2550 -lt 10000 ']'
+ FICHIER=velocity_Z_it002550
+ '[' 2550 -lt 1000 ']'
+ '[' 2550 -lt 100 ']'
+ echo velocity_Z_it002550
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002550 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2550+25
+ '[' 2575 -le 8500 ']'
+ '[' 2575 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2575
+ '[' 2575 -lt 100000 ']'
+ FICHIER=velocity_Z_it02575
+ '[' 2575 -lt 10000 ']'
+ FICHIER=velocity_Z_it002575
+ '[' 2575 -lt 1000 ']'
+ '[' 2575 -lt 100 ']'
+ echo velocity_Z_it002575
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002575 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2575+25
+ '[' 2600 -le 8500 ']'
+ '[' 2600 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2600
+ '[' 2600 -lt 100000 ']'
+ FICHIER=velocity_Z_it02600
+ '[' 2600 -lt 10000 ']'
+ FICHIER=velocity_Z_it002600
+ '[' 2600 -lt 1000 ']'
+ '[' 2600 -lt 100 ']'
+ echo velocity_Z_it002600
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002600 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2600+25
+ '[' 2625 -le 8500 ']'
+ '[' 2625 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2625
+ '[' 2625 -lt 100000 ']'
+ FICHIER=velocity_Z_it02625
+ '[' 2625 -lt 10000 ']'
+ FICHIER=velocity_Z_it002625
+ '[' 2625 -lt 1000 ']'
+ '[' 2625 -lt 100 ']'
+ echo velocity_Z_it002625
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002625 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2625+25
+ '[' 2650 -le 8500 ']'
+ '[' 2650 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2650
+ '[' 2650 -lt 100000 ']'
+ FICHIER=velocity_Z_it02650
+ '[' 2650 -lt 10000 ']'
+ FICHIER=velocity_Z_it002650
+ '[' 2650 -lt 1000 ']'
+ '[' 2650 -lt 100 ']'
+ echo velocity_Z_it002650
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002650 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2650+25
+ '[' 2675 -le 8500 ']'
+ '[' 2675 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2675
+ '[' 2675 -lt 100000 ']'
+ FICHIER=velocity_Z_it02675
+ '[' 2675 -lt 10000 ']'
+ FICHIER=velocity_Z_it002675
+ '[' 2675 -lt 1000 ']'
+ '[' 2675 -lt 100 ']'
+ echo velocity_Z_it002675
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002675 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2675+25
+ '[' 2700 -le 8500 ']'
+ '[' 2700 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2700
+ '[' 2700 -lt 100000 ']'
+ FICHIER=velocity_Z_it02700
+ '[' 2700 -lt 10000 ']'
+ FICHIER=velocity_Z_it002700
+ '[' 2700 -lt 1000 ']'
+ '[' 2700 -lt 100 ']'
+ echo velocity_Z_it002700
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002700 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2700+25
+ '[' 2725 -le 8500 ']'
+ '[' 2725 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2725
+ '[' 2725 -lt 100000 ']'
+ FICHIER=velocity_Z_it02725
+ '[' 2725 -lt 10000 ']'
+ FICHIER=velocity_Z_it002725
+ '[' 2725 -lt 1000 ']'
+ '[' 2725 -lt 100 ']'
+ echo velocity_Z_it002725
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002725 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2725+25
+ '[' 2750 -le 8500 ']'
+ '[' 2750 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2750
+ '[' 2750 -lt 100000 ']'
+ FICHIER=velocity_Z_it02750
+ '[' 2750 -lt 10000 ']'
+ FICHIER=velocity_Z_it002750
+ '[' 2750 -lt 1000 ']'
+ '[' 2750 -lt 100 ']'
+ echo velocity_Z_it002750
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002750 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2750+25
+ '[' 2775 -le 8500 ']'
+ '[' 2775 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2775
+ '[' 2775 -lt 100000 ']'
+ FICHIER=velocity_Z_it02775
+ '[' 2775 -lt 10000 ']'
+ FICHIER=velocity_Z_it002775
+ '[' 2775 -lt 1000 ']'
+ '[' 2775 -lt 100 ']'
+ echo velocity_Z_it002775
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002775 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2775+25
+ '[' 2800 -le 8500 ']'
+ '[' 2800 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2800
+ '[' 2800 -lt 100000 ']'
+ FICHIER=velocity_Z_it02800
+ '[' 2800 -lt 10000 ']'
+ FICHIER=velocity_Z_it002800
+ '[' 2800 -lt 1000 ']'
+ '[' 2800 -lt 100 ']'
+ echo velocity_Z_it002800
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002800 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2800+25
+ '[' 2825 -le 8500 ']'
+ '[' 2825 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2825
+ '[' 2825 -lt 100000 ']'
+ FICHIER=velocity_Z_it02825
+ '[' 2825 -lt 10000 ']'
+ FICHIER=velocity_Z_it002825
+ '[' 2825 -lt 1000 ']'
+ '[' 2825 -lt 100 ']'
+ echo velocity_Z_it002825
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002825 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2825+25
+ '[' 2850 -le 8500 ']'
+ '[' 2850 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2850
+ '[' 2850 -lt 100000 ']'
+ FICHIER=velocity_Z_it02850
+ '[' 2850 -lt 10000 ']'
+ FICHIER=velocity_Z_it002850
+ '[' 2850 -lt 1000 ']'
+ '[' 2850 -lt 100 ']'
+ echo velocity_Z_it002850
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002850 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2850+25
+ '[' 2875 -le 8500 ']'
+ '[' 2875 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2875
+ '[' 2875 -lt 100000 ']'
+ FICHIER=velocity_Z_it02875
+ '[' 2875 -lt 10000 ']'
+ FICHIER=velocity_Z_it002875
+ '[' 2875 -lt 1000 ']'
+ '[' 2875 -lt 100 ']'
+ echo velocity_Z_it002875
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002875 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2875+25
+ '[' 2900 -le 8500 ']'
+ '[' 2900 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2900
+ '[' 2900 -lt 100000 ']'
+ FICHIER=velocity_Z_it02900
+ '[' 2900 -lt 10000 ']'
+ FICHIER=velocity_Z_it002900
+ '[' 2900 -lt 1000 ']'
+ '[' 2900 -lt 100 ']'
+ echo velocity_Z_it002900
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002900 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2900+25
+ '[' 2925 -le 8500 ']'
+ '[' 2925 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2925
+ '[' 2925 -lt 100000 ']'
+ FICHIER=velocity_Z_it02925
+ '[' 2925 -lt 10000 ']'
+ FICHIER=velocity_Z_it002925
+ '[' 2925 -lt 1000 ']'
+ '[' 2925 -lt 100 ']'
+ echo velocity_Z_it002925
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002925 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2925+25
+ '[' 2950 -le 8500 ']'
+ '[' 2950 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2950
+ '[' 2950 -lt 100000 ']'
+ FICHIER=velocity_Z_it02950
+ '[' 2950 -lt 10000 ']'
+ FICHIER=velocity_Z_it002950
+ '[' 2950 -lt 1000 ']'
+ '[' 2950 -lt 100 ']'
+ echo velocity_Z_it002950
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002950 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2950+25
+ '[' 2975 -le 8500 ']'
+ '[' 2975 -lt 1000000 ']'
+ FICHIER=velocity_Z_it2975
+ '[' 2975 -lt 100000 ']'
+ FICHIER=velocity_Z_it02975
+ '[' 2975 -lt 10000 ']'
+ FICHIER=velocity_Z_it002975
+ '[' 2975 -lt 1000 ']'
+ '[' 2975 -lt 100 ']'
+ echo velocity_Z_it002975
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it002975 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=2975+25
+ '[' 3000 -le 8500 ']'
+ '[' 3000 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3000
+ '[' 3000 -lt 100000 ']'
+ FICHIER=velocity_Z_it03000
+ '[' 3000 -lt 10000 ']'
+ FICHIER=velocity_Z_it003000
+ '[' 3000 -lt 1000 ']'
+ '[' 3000 -lt 100 ']'
+ echo velocity_Z_it003000
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003000 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3000+25
+ '[' 3025 -le 8500 ']'
+ '[' 3025 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3025
+ '[' 3025 -lt 100000 ']'
+ FICHIER=velocity_Z_it03025
+ '[' 3025 -lt 10000 ']'
+ FICHIER=velocity_Z_it003025
+ '[' 3025 -lt 1000 ']'
+ '[' 3025 -lt 100 ']'
+ echo velocity_Z_it003025
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003025 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3025+25
+ '[' 3050 -le 8500 ']'
+ '[' 3050 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3050
+ '[' 3050 -lt 100000 ']'
+ FICHIER=velocity_Z_it03050
+ '[' 3050 -lt 10000 ']'
+ FICHIER=velocity_Z_it003050
+ '[' 3050 -lt 1000 ']'
+ '[' 3050 -lt 100 ']'
+ echo velocity_Z_it003050
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003050 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3050+25
+ '[' 3075 -le 8500 ']'
+ '[' 3075 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3075
+ '[' 3075 -lt 100000 ']'
+ FICHIER=velocity_Z_it03075
+ '[' 3075 -lt 10000 ']'
+ FICHIER=velocity_Z_it003075
+ '[' 3075 -lt 1000 ']'
+ '[' 3075 -lt 100 ']'
+ echo velocity_Z_it003075
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003075 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3075+25
+ '[' 3100 -le 8500 ']'
+ '[' 3100 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3100
+ '[' 3100 -lt 100000 ']'
+ FICHIER=velocity_Z_it03100
+ '[' 3100 -lt 10000 ']'
+ FICHIER=velocity_Z_it003100
+ '[' 3100 -lt 1000 ']'
+ '[' 3100 -lt 100 ']'
+ echo velocity_Z_it003100
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003100 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3100+25
+ '[' 3125 -le 8500 ']'
+ '[' 3125 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3125
+ '[' 3125 -lt 100000 ']'
+ FICHIER=velocity_Z_it03125
+ '[' 3125 -lt 10000 ']'
+ FICHIER=velocity_Z_it003125
+ '[' 3125 -lt 1000 ']'
+ '[' 3125 -lt 100 ']'
+ echo velocity_Z_it003125
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003125 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3125+25
+ '[' 3150 -le 8500 ']'
+ '[' 3150 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3150
+ '[' 3150 -lt 100000 ']'
+ FICHIER=velocity_Z_it03150
+ '[' 3150 -lt 10000 ']'
+ FICHIER=velocity_Z_it003150
+ '[' 3150 -lt 1000 ']'
+ '[' 3150 -lt 100 ']'
+ echo velocity_Z_it003150
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003150 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3150+25
+ '[' 3175 -le 8500 ']'
+ '[' 3175 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3175
+ '[' 3175 -lt 100000 ']'
+ FICHIER=velocity_Z_it03175
+ '[' 3175 -lt 10000 ']'
+ FICHIER=velocity_Z_it003175
+ '[' 3175 -lt 1000 ']'
+ '[' 3175 -lt 100 ']'
+ echo velocity_Z_it003175
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003175 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3175+25
+ '[' 3200 -le 8500 ']'
+ '[' 3200 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3200
+ '[' 3200 -lt 100000 ']'
+ FICHIER=velocity_Z_it03200
+ '[' 3200 -lt 10000 ']'
+ FICHIER=velocity_Z_it003200
+ '[' 3200 -lt 1000 ']'
+ '[' 3200 -lt 100 ']'
+ echo velocity_Z_it003200
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003200 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3200+25
+ '[' 3225 -le 8500 ']'
+ '[' 3225 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3225
+ '[' 3225 -lt 100000 ']'
+ FICHIER=velocity_Z_it03225
+ '[' 3225 -lt 10000 ']'
+ FICHIER=velocity_Z_it003225
+ '[' 3225 -lt 1000 ']'
+ '[' 3225 -lt 100 ']'
+ echo velocity_Z_it003225
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003225 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3225+25
+ '[' 3250 -le 8500 ']'
+ '[' 3250 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3250
+ '[' 3250 -lt 100000 ']'
+ FICHIER=velocity_Z_it03250
+ '[' 3250 -lt 10000 ']'
+ FICHIER=velocity_Z_it003250
+ '[' 3250 -lt 1000 ']'
+ '[' 3250 -lt 100 ']'
+ echo velocity_Z_it003250
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003250 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3250+25
+ '[' 3275 -le 8500 ']'
+ '[' 3275 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3275
+ '[' 3275 -lt 100000 ']'
+ FICHIER=velocity_Z_it03275
+ '[' 3275 -lt 10000 ']'
+ FICHIER=velocity_Z_it003275
+ '[' 3275 -lt 1000 ']'
+ '[' 3275 -lt 100 ']'
+ echo velocity_Z_it003275
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003275 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3275+25
+ '[' 3300 -le 8500 ']'
+ '[' 3300 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3300
+ '[' 3300 -lt 100000 ']'
+ FICHIER=velocity_Z_it03300
+ '[' 3300 -lt 10000 ']'
+ FICHIER=velocity_Z_it003300
+ '[' 3300 -lt 1000 ']'
+ '[' 3300 -lt 100 ']'
+ echo velocity_Z_it003300
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003300 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3300+25
+ '[' 3325 -le 8500 ']'
+ '[' 3325 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3325
+ '[' 3325 -lt 100000 ']'
+ FICHIER=velocity_Z_it03325
+ '[' 3325 -lt 10000 ']'
+ FICHIER=velocity_Z_it003325
+ '[' 3325 -lt 1000 ']'
+ '[' 3325 -lt 100 ']'
+ echo velocity_Z_it003325
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003325 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3325+25
+ '[' 3350 -le 8500 ']'
+ '[' 3350 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3350
+ '[' 3350 -lt 100000 ']'
+ FICHIER=velocity_Z_it03350
+ '[' 3350 -lt 10000 ']'
+ FICHIER=velocity_Z_it003350
+ '[' 3350 -lt 1000 ']'
+ '[' 3350 -lt 100 ']'
+ echo velocity_Z_it003350
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003350 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3350+25
+ '[' 3375 -le 8500 ']'
+ '[' 3375 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3375
+ '[' 3375 -lt 100000 ']'
+ FICHIER=velocity_Z_it03375
+ '[' 3375 -lt 10000 ']'
+ FICHIER=velocity_Z_it003375
+ '[' 3375 -lt 1000 ']'
+ '[' 3375 -lt 100 ']'
+ echo velocity_Z_it003375
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003375 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3375+25
+ '[' 3400 -le 8500 ']'
+ '[' 3400 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3400
+ '[' 3400 -lt 100000 ']'
+ FICHIER=velocity_Z_it03400
+ '[' 3400 -lt 10000 ']'
+ FICHIER=velocity_Z_it003400
+ '[' 3400 -lt 1000 ']'
+ '[' 3400 -lt 100 ']'
+ echo velocity_Z_it003400
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003400 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3400+25
+ '[' 3425 -le 8500 ']'
+ '[' 3425 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3425
+ '[' 3425 -lt 100000 ']'
+ FICHIER=velocity_Z_it03425
+ '[' 3425 -lt 10000 ']'
+ FICHIER=velocity_Z_it003425
+ '[' 3425 -lt 1000 ']'
+ '[' 3425 -lt 100 ']'
+ echo velocity_Z_it003425
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003425 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3425+25
+ '[' 3450 -le 8500 ']'
+ '[' 3450 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3450
+ '[' 3450 -lt 100000 ']'
+ FICHIER=velocity_Z_it03450
+ '[' 3450 -lt 10000 ']'
+ FICHIER=velocity_Z_it003450
+ '[' 3450 -lt 1000 ']'
+ '[' 3450 -lt 100 ']'
+ echo velocity_Z_it003450
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003450 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3450+25
+ '[' 3475 -le 8500 ']'
+ '[' 3475 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3475
+ '[' 3475 -lt 100000 ']'
+ FICHIER=velocity_Z_it03475
+ '[' 3475 -lt 10000 ']'
+ FICHIER=velocity_Z_it003475
+ '[' 3475 -lt 1000 ']'
+ '[' 3475 -lt 100 ']'
+ echo velocity_Z_it003475
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003475 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3475+25
+ '[' 3500 -le 8500 ']'
+ '[' 3500 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3500
+ '[' 3500 -lt 100000 ']'
+ FICHIER=velocity_Z_it03500
+ '[' 3500 -lt 10000 ']'
+ FICHIER=velocity_Z_it003500
+ '[' 3500 -lt 1000 ']'
+ '[' 3500 -lt 100 ']'
+ echo velocity_Z_it003500
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003500 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3500+25
+ '[' 3525 -le 8500 ']'
+ '[' 3525 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3525
+ '[' 3525 -lt 100000 ']'
+ FICHIER=velocity_Z_it03525
+ '[' 3525 -lt 10000 ']'
+ FICHIER=velocity_Z_it003525
+ '[' 3525 -lt 1000 ']'
+ '[' 3525 -lt 100 ']'
+ echo velocity_Z_it003525
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003525 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3525+25
+ '[' 3550 -le 8500 ']'
+ '[' 3550 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3550
+ '[' 3550 -lt 100000 ']'
+ FICHIER=velocity_Z_it03550
+ '[' 3550 -lt 10000 ']'
+ FICHIER=velocity_Z_it003550
+ '[' 3550 -lt 1000 ']'
+ '[' 3550 -lt 100 ']'
+ echo velocity_Z_it003550
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003550 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3550+25
+ '[' 3575 -le 8500 ']'
+ '[' 3575 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3575
+ '[' 3575 -lt 100000 ']'
+ FICHIER=velocity_Z_it03575
+ '[' 3575 -lt 10000 ']'
+ FICHIER=velocity_Z_it003575
+ '[' 3575 -lt 1000 ']'
+ '[' 3575 -lt 100 ']'
+ echo velocity_Z_it003575
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003575 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3575+25
+ '[' 3600 -le 8500 ']'
+ '[' 3600 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3600
+ '[' 3600 -lt 100000 ']'
+ FICHIER=velocity_Z_it03600
+ '[' 3600 -lt 10000 ']'
+ FICHIER=velocity_Z_it003600
+ '[' 3600 -lt 1000 ']'
+ '[' 3600 -lt 100 ']'
+ echo velocity_Z_it003600
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003600 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3600+25
+ '[' 3625 -le 8500 ']'
+ '[' 3625 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3625
+ '[' 3625 -lt 100000 ']'
+ FICHIER=velocity_Z_it03625
+ '[' 3625 -lt 10000 ']'
+ FICHIER=velocity_Z_it003625
+ '[' 3625 -lt 1000 ']'
+ '[' 3625 -lt 100 ']'
+ echo velocity_Z_it003625
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003625 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3625+25
+ '[' 3650 -le 8500 ']'
+ '[' 3650 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3650
+ '[' 3650 -lt 100000 ']'
+ FICHIER=velocity_Z_it03650
+ '[' 3650 -lt 10000 ']'
+ FICHIER=velocity_Z_it003650
+ '[' 3650 -lt 1000 ']'
+ '[' 3650 -lt 100 ']'
+ echo velocity_Z_it003650
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003650 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3650+25
+ '[' 3675 -le 8500 ']'
+ '[' 3675 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3675
+ '[' 3675 -lt 100000 ']'
+ FICHIER=velocity_Z_it03675
+ '[' 3675 -lt 10000 ']'
+ FICHIER=velocity_Z_it003675
+ '[' 3675 -lt 1000 ']'
+ '[' 3675 -lt 100 ']'
+ echo velocity_Z_it003675
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003675 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3675+25
+ '[' 3700 -le 8500 ']'
+ '[' 3700 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3700
+ '[' 3700 -lt 100000 ']'
+ FICHIER=velocity_Z_it03700
+ '[' 3700 -lt 10000 ']'
+ FICHIER=velocity_Z_it003700
+ '[' 3700 -lt 1000 ']'
+ '[' 3700 -lt 100 ']'
+ echo velocity_Z_it003700
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003700 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3700+25
+ '[' 3725 -le 8500 ']'
+ '[' 3725 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3725
+ '[' 3725 -lt 100000 ']'
+ FICHIER=velocity_Z_it03725
+ '[' 3725 -lt 10000 ']'
+ FICHIER=velocity_Z_it003725
+ '[' 3725 -lt 1000 ']'
+ '[' 3725 -lt 100 ']'
+ echo velocity_Z_it003725
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003725 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3725+25
+ '[' 3750 -le 8500 ']'
+ '[' 3750 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3750
+ '[' 3750 -lt 100000 ']'
+ FICHIER=velocity_Z_it03750
+ '[' 3750 -lt 10000 ']'
+ FICHIER=velocity_Z_it003750
+ '[' 3750 -lt 1000 ']'
+ '[' 3750 -lt 100 ']'
+ echo velocity_Z_it003750
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003750 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3750+25
+ '[' 3775 -le 8500 ']'
+ '[' 3775 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3775
+ '[' 3775 -lt 100000 ']'
+ FICHIER=velocity_Z_it03775
+ '[' 3775 -lt 10000 ']'
+ FICHIER=velocity_Z_it003775
+ '[' 3775 -lt 1000 ']'
+ '[' 3775 -lt 100 ']'
+ echo velocity_Z_it003775
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003775 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3775+25
+ '[' 3800 -le 8500 ']'
+ '[' 3800 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3800
+ '[' 3800 -lt 100000 ']'
+ FICHIER=velocity_Z_it03800
+ '[' 3800 -lt 10000 ']'
+ FICHIER=velocity_Z_it003800
+ '[' 3800 -lt 1000 ']'
+ '[' 3800 -lt 100 ']'
+ echo velocity_Z_it003800
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003800 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3800+25
+ '[' 3825 -le 8500 ']'
+ '[' 3825 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3825
+ '[' 3825 -lt 100000 ']'
+ FICHIER=velocity_Z_it03825
+ '[' 3825 -lt 10000 ']'
+ FICHIER=velocity_Z_it003825
+ '[' 3825 -lt 1000 ']'
+ '[' 3825 -lt 100 ']'
+ echo velocity_Z_it003825
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003825 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3825+25
+ '[' 3850 -le 8500 ']'
+ '[' 3850 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3850
+ '[' 3850 -lt 100000 ']'
+ FICHIER=velocity_Z_it03850
+ '[' 3850 -lt 10000 ']'
+ FICHIER=velocity_Z_it003850
+ '[' 3850 -lt 1000 ']'
+ '[' 3850 -lt 100 ']'
+ echo velocity_Z_it003850
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003850 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3850+25
+ '[' 3875 -le 8500 ']'
+ '[' 3875 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3875
+ '[' 3875 -lt 100000 ']'
+ FICHIER=velocity_Z_it03875
+ '[' 3875 -lt 10000 ']'
+ FICHIER=velocity_Z_it003875
+ '[' 3875 -lt 1000 ']'
+ '[' 3875 -lt 100 ']'
+ echo velocity_Z_it003875
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003875 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3875+25
+ '[' 3900 -le 8500 ']'
+ '[' 3900 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3900
+ '[' 3900 -lt 100000 ']'
+ FICHIER=velocity_Z_it03900
+ '[' 3900 -lt 10000 ']'
+ FICHIER=velocity_Z_it003900
+ '[' 3900 -lt 1000 ']'
+ '[' 3900 -lt 100 ']'
+ echo velocity_Z_it003900
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003900 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3900+25
+ '[' 3925 -le 8500 ']'
+ '[' 3925 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3925
+ '[' 3925 -lt 100000 ']'
+ FICHIER=velocity_Z_it03925
+ '[' 3925 -lt 10000 ']'
+ FICHIER=velocity_Z_it003925
+ '[' 3925 -lt 1000 ']'
+ '[' 3925 -lt 100 ']'
+ echo velocity_Z_it003925
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003925 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3925+25
+ '[' 3950 -le 8500 ']'
+ '[' 3950 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3950
+ '[' 3950 -lt 100000 ']'
+ FICHIER=velocity_Z_it03950
+ '[' 3950 -lt 10000 ']'
+ FICHIER=velocity_Z_it003950
+ '[' 3950 -lt 1000 ']'
+ '[' 3950 -lt 100 ']'
+ echo velocity_Z_it003950
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003950 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3950+25
+ '[' 3975 -le 8500 ']'
+ '[' 3975 -lt 1000000 ']'
+ FICHIER=velocity_Z_it3975
+ '[' 3975 -lt 100000 ']'
+ FICHIER=velocity_Z_it03975
+ '[' 3975 -lt 10000 ']'
+ FICHIER=velocity_Z_it003975
+ '[' 3975 -lt 1000 ']'
+ '[' 3975 -lt 100 ']'
+ echo velocity_Z_it003975
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it003975 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=3975+25
+ '[' 4000 -le 8500 ']'
+ '[' 4000 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4000
+ '[' 4000 -lt 100000 ']'
+ FICHIER=velocity_Z_it04000
+ '[' 4000 -lt 10000 ']'
+ FICHIER=velocity_Z_it004000
+ '[' 4000 -lt 1000 ']'
+ '[' 4000 -lt 100 ']'
+ echo velocity_Z_it004000
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004000 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4000+25
+ '[' 4025 -le 8500 ']'
+ '[' 4025 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4025
+ '[' 4025 -lt 100000 ']'
+ FICHIER=velocity_Z_it04025
+ '[' 4025 -lt 10000 ']'
+ FICHIER=velocity_Z_it004025
+ '[' 4025 -lt 1000 ']'
+ '[' 4025 -lt 100 ']'
+ echo velocity_Z_it004025
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004025 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4025+25
+ '[' 4050 -le 8500 ']'
+ '[' 4050 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4050
+ '[' 4050 -lt 100000 ']'
+ FICHIER=velocity_Z_it04050
+ '[' 4050 -lt 10000 ']'
+ FICHIER=velocity_Z_it004050
+ '[' 4050 -lt 1000 ']'
+ '[' 4050 -lt 100 ']'
+ echo velocity_Z_it004050
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004050 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4050+25
+ '[' 4075 -le 8500 ']'
+ '[' 4075 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4075
+ '[' 4075 -lt 100000 ']'
+ FICHIER=velocity_Z_it04075
+ '[' 4075 -lt 10000 ']'
+ FICHIER=velocity_Z_it004075
+ '[' 4075 -lt 1000 ']'
+ '[' 4075 -lt 100 ']'
+ echo velocity_Z_it004075
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004075 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4075+25
+ '[' 4100 -le 8500 ']'
+ '[' 4100 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4100
+ '[' 4100 -lt 100000 ']'
+ FICHIER=velocity_Z_it04100
+ '[' 4100 -lt 10000 ']'
+ FICHIER=velocity_Z_it004100
+ '[' 4100 -lt 1000 ']'
+ '[' 4100 -lt 100 ']'
+ echo velocity_Z_it004100
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004100 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4100+25
+ '[' 4125 -le 8500 ']'
+ '[' 4125 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4125
+ '[' 4125 -lt 100000 ']'
+ FICHIER=velocity_Z_it04125
+ '[' 4125 -lt 10000 ']'
+ FICHIER=velocity_Z_it004125
+ '[' 4125 -lt 1000 ']'
+ '[' 4125 -lt 100 ']'
+ echo velocity_Z_it004125
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004125 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4125+25
+ '[' 4150 -le 8500 ']'
+ '[' 4150 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4150
+ '[' 4150 -lt 100000 ']'
+ FICHIER=velocity_Z_it04150
+ '[' 4150 -lt 10000 ']'
+ FICHIER=velocity_Z_it004150
+ '[' 4150 -lt 1000 ']'
+ '[' 4150 -lt 100 ']'
+ echo velocity_Z_it004150
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004150 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4150+25
+ '[' 4175 -le 8500 ']'
+ '[' 4175 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4175
+ '[' 4175 -lt 100000 ']'
+ FICHIER=velocity_Z_it04175
+ '[' 4175 -lt 10000 ']'
+ FICHIER=velocity_Z_it004175
+ '[' 4175 -lt 1000 ']'
+ '[' 4175 -lt 100 ']'
+ echo velocity_Z_it004175
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004175 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4175+25
+ '[' 4200 -le 8500 ']'
+ '[' 4200 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4200
+ '[' 4200 -lt 100000 ']'
+ FICHIER=velocity_Z_it04200
+ '[' 4200 -lt 10000 ']'
+ FICHIER=velocity_Z_it004200
+ '[' 4200 -lt 1000 ']'
+ '[' 4200 -lt 100 ']'
+ echo velocity_Z_it004200
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004200 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4200+25
+ '[' 4225 -le 8500 ']'
+ '[' 4225 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4225
+ '[' 4225 -lt 100000 ']'
+ FICHIER=velocity_Z_it04225
+ '[' 4225 -lt 10000 ']'
+ FICHIER=velocity_Z_it004225
+ '[' 4225 -lt 1000 ']'
+ '[' 4225 -lt 100 ']'
+ echo velocity_Z_it004225
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004225 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4225+25
+ '[' 4250 -le 8500 ']'
+ '[' 4250 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4250
+ '[' 4250 -lt 100000 ']'
+ FICHIER=velocity_Z_it04250
+ '[' 4250 -lt 10000 ']'
+ FICHIER=velocity_Z_it004250
+ '[' 4250 -lt 1000 ']'
+ '[' 4250 -lt 100 ']'
+ echo velocity_Z_it004250
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004250 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4250+25
+ '[' 4275 -le 8500 ']'
+ '[' 4275 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4275
+ '[' 4275 -lt 100000 ']'
+ FICHIER=velocity_Z_it04275
+ '[' 4275 -lt 10000 ']'
+ FICHIER=velocity_Z_it004275
+ '[' 4275 -lt 1000 ']'
+ '[' 4275 -lt 100 ']'
+ echo velocity_Z_it004275
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004275 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4275+25
+ '[' 4300 -le 8500 ']'
+ '[' 4300 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4300
+ '[' 4300 -lt 100000 ']'
+ FICHIER=velocity_Z_it04300
+ '[' 4300 -lt 10000 ']'
+ FICHIER=velocity_Z_it004300
+ '[' 4300 -lt 1000 ']'
+ '[' 4300 -lt 100 ']'
+ echo velocity_Z_it004300
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004300 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4300+25
+ '[' 4325 -le 8500 ']'
+ '[' 4325 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4325
+ '[' 4325 -lt 100000 ']'
+ FICHIER=velocity_Z_it04325
+ '[' 4325 -lt 10000 ']'
+ FICHIER=velocity_Z_it004325
+ '[' 4325 -lt 1000 ']'
+ '[' 4325 -lt 100 ']'
+ echo velocity_Z_it004325
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004325 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4325+25
+ '[' 4350 -le 8500 ']'
+ '[' 4350 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4350
+ '[' 4350 -lt 100000 ']'
+ FICHIER=velocity_Z_it04350
+ '[' 4350 -lt 10000 ']'
+ FICHIER=velocity_Z_it004350
+ '[' 4350 -lt 1000 ']'
+ '[' 4350 -lt 100 ']'
+ echo velocity_Z_it004350
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004350 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4350+25
+ '[' 4375 -le 8500 ']'
+ '[' 4375 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4375
+ '[' 4375 -lt 100000 ']'
+ FICHIER=velocity_Z_it04375
+ '[' 4375 -lt 10000 ']'
+ FICHIER=velocity_Z_it004375
+ '[' 4375 -lt 1000 ']'
+ '[' 4375 -lt 100 ']'
+ echo velocity_Z_it004375
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004375 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4375+25
+ '[' 4400 -le 8500 ']'
+ '[' 4400 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4400
+ '[' 4400 -lt 100000 ']'
+ FICHIER=velocity_Z_it04400
+ '[' 4400 -lt 10000 ']'
+ FICHIER=velocity_Z_it004400
+ '[' 4400 -lt 1000 ']'
+ '[' 4400 -lt 100 ']'
+ echo velocity_Z_it004400
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004400 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4400+25
+ '[' 4425 -le 8500 ']'
+ '[' 4425 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4425
+ '[' 4425 -lt 100000 ']'
+ FICHIER=velocity_Z_it04425
+ '[' 4425 -lt 10000 ']'
+ FICHIER=velocity_Z_it004425
+ '[' 4425 -lt 1000 ']'
+ '[' 4425 -lt 100 ']'
+ echo velocity_Z_it004425
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004425 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4425+25
+ '[' 4450 -le 8500 ']'
+ '[' 4450 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4450
+ '[' 4450 -lt 100000 ']'
+ FICHIER=velocity_Z_it04450
+ '[' 4450 -lt 10000 ']'
+ FICHIER=velocity_Z_it004450
+ '[' 4450 -lt 1000 ']'
+ '[' 4450 -lt 100 ']'
+ echo velocity_Z_it004450
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004450 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4450+25
+ '[' 4475 -le 8500 ']'
+ '[' 4475 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4475
+ '[' 4475 -lt 100000 ']'
+ FICHIER=velocity_Z_it04475
+ '[' 4475 -lt 10000 ']'
+ FICHIER=velocity_Z_it004475
+ '[' 4475 -lt 1000 ']'
+ '[' 4475 -lt 100 ']'
+ echo velocity_Z_it004475
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004475 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4475+25
+ '[' 4500 -le 8500 ']'
+ '[' 4500 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4500
+ '[' 4500 -lt 100000 ']'
+ FICHIER=velocity_Z_it04500
+ '[' 4500 -lt 10000 ']'
+ FICHIER=velocity_Z_it004500
+ '[' 4500 -lt 1000 ']'
+ '[' 4500 -lt 100 ']'
+ echo velocity_Z_it004500
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004500 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4500+25
+ '[' 4525 -le 8500 ']'
+ '[' 4525 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4525
+ '[' 4525 -lt 100000 ']'
+ FICHIER=velocity_Z_it04525
+ '[' 4525 -lt 10000 ']'
+ FICHIER=velocity_Z_it004525
+ '[' 4525 -lt 1000 ']'
+ '[' 4525 -lt 100 ']'
+ echo velocity_Z_it004525
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004525 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4525+25
+ '[' 4550 -le 8500 ']'
+ '[' 4550 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4550
+ '[' 4550 -lt 100000 ']'
+ FICHIER=velocity_Z_it04550
+ '[' 4550 -lt 10000 ']'
+ FICHIER=velocity_Z_it004550
+ '[' 4550 -lt 1000 ']'
+ '[' 4550 -lt 100 ']'
+ echo velocity_Z_it004550
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004550 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4550+25
+ '[' 4575 -le 8500 ']'
+ '[' 4575 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4575
+ '[' 4575 -lt 100000 ']'
+ FICHIER=velocity_Z_it04575
+ '[' 4575 -lt 10000 ']'
+ FICHIER=velocity_Z_it004575
+ '[' 4575 -lt 1000 ']'
+ '[' 4575 -lt 100 ']'
+ echo velocity_Z_it004575
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004575 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4575+25
+ '[' 4600 -le 8500 ']'
+ '[' 4600 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4600
+ '[' 4600 -lt 100000 ']'
+ FICHIER=velocity_Z_it04600
+ '[' 4600 -lt 10000 ']'
+ FICHIER=velocity_Z_it004600
+ '[' 4600 -lt 1000 ']'
+ '[' 4600 -lt 100 ']'
+ echo velocity_Z_it004600
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004600 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4600+25
+ '[' 4625 -le 8500 ']'
+ '[' 4625 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4625
+ '[' 4625 -lt 100000 ']'
+ FICHIER=velocity_Z_it04625
+ '[' 4625 -lt 10000 ']'
+ FICHIER=velocity_Z_it004625
+ '[' 4625 -lt 1000 ']'
+ '[' 4625 -lt 100 ']'
+ echo velocity_Z_it004625
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004625 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4625+25
+ '[' 4650 -le 8500 ']'
+ '[' 4650 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4650
+ '[' 4650 -lt 100000 ']'
+ FICHIER=velocity_Z_it04650
+ '[' 4650 -lt 10000 ']'
+ FICHIER=velocity_Z_it004650
+ '[' 4650 -lt 1000 ']'
+ '[' 4650 -lt 100 ']'
+ echo velocity_Z_it004650
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004650 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4650+25
+ '[' 4675 -le 8500 ']'
+ '[' 4675 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4675
+ '[' 4675 -lt 100000 ']'
+ FICHIER=velocity_Z_it04675
+ '[' 4675 -lt 10000 ']'
+ FICHIER=velocity_Z_it004675
+ '[' 4675 -lt 1000 ']'
+ '[' 4675 -lt 100 ']'
+ echo velocity_Z_it004675
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004675 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4675+25
+ '[' 4700 -le 8500 ']'
+ '[' 4700 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4700
+ '[' 4700 -lt 100000 ']'
+ FICHIER=velocity_Z_it04700
+ '[' 4700 -lt 10000 ']'
+ FICHIER=velocity_Z_it004700
+ '[' 4700 -lt 1000 ']'
+ '[' 4700 -lt 100 ']'
+ echo velocity_Z_it004700
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004700 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4700+25
+ '[' 4725 -le 8500 ']'
+ '[' 4725 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4725
+ '[' 4725 -lt 100000 ']'
+ FICHIER=velocity_Z_it04725
+ '[' 4725 -lt 10000 ']'
+ FICHIER=velocity_Z_it004725
+ '[' 4725 -lt 1000 ']'
+ '[' 4725 -lt 100 ']'
+ echo velocity_Z_it004725
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004725 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4725+25
+ '[' 4750 -le 8500 ']'
+ '[' 4750 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4750
+ '[' 4750 -lt 100000 ']'
+ FICHIER=velocity_Z_it04750
+ '[' 4750 -lt 10000 ']'
+ FICHIER=velocity_Z_it004750
+ '[' 4750 -lt 1000 ']'
+ '[' 4750 -lt 100 ']'
+ echo velocity_Z_it004750
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004750 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4750+25
+ '[' 4775 -le 8500 ']'
+ '[' 4775 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4775
+ '[' 4775 -lt 100000 ']'
+ FICHIER=velocity_Z_it04775
+ '[' 4775 -lt 10000 ']'
+ FICHIER=velocity_Z_it004775
+ '[' 4775 -lt 1000 ']'
+ '[' 4775 -lt 100 ']'
+ echo velocity_Z_it004775
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004775 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4775+25
+ '[' 4800 -le 8500 ']'
+ '[' 4800 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4800
+ '[' 4800 -lt 100000 ']'
+ FICHIER=velocity_Z_it04800
+ '[' 4800 -lt 10000 ']'
+ FICHIER=velocity_Z_it004800
+ '[' 4800 -lt 1000 ']'
+ '[' 4800 -lt 100 ']'
+ echo velocity_Z_it004800
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004800 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4800+25
+ '[' 4825 -le 8500 ']'
+ '[' 4825 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4825
+ '[' 4825 -lt 100000 ']'
+ FICHIER=velocity_Z_it04825
+ '[' 4825 -lt 10000 ']'
+ FICHIER=velocity_Z_it004825
+ '[' 4825 -lt 1000 ']'
+ '[' 4825 -lt 100 ']'
+ echo velocity_Z_it004825
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004825 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4825+25
+ '[' 4850 -le 8500 ']'
+ '[' 4850 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4850
+ '[' 4850 -lt 100000 ']'
+ FICHIER=velocity_Z_it04850
+ '[' 4850 -lt 10000 ']'
+ FICHIER=velocity_Z_it004850
+ '[' 4850 -lt 1000 ']'
+ '[' 4850 -lt 100 ']'
+ echo velocity_Z_it004850
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004850 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4850+25
+ '[' 4875 -le 8500 ']'
+ '[' 4875 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4875
+ '[' 4875 -lt 100000 ']'
+ FICHIER=velocity_Z_it04875
+ '[' 4875 -lt 10000 ']'
+ FICHIER=velocity_Z_it004875
+ '[' 4875 -lt 1000 ']'
+ '[' 4875 -lt 100 ']'
+ echo velocity_Z_it004875
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004875 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4875+25
+ '[' 4900 -le 8500 ']'
+ '[' 4900 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4900
+ '[' 4900 -lt 100000 ']'
+ FICHIER=velocity_Z_it04900
+ '[' 4900 -lt 10000 ']'
+ FICHIER=velocity_Z_it004900
+ '[' 4900 -lt 1000 ']'
+ '[' 4900 -lt 100 ']'
+ echo velocity_Z_it004900
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004900 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4900+25
+ '[' 4925 -le 8500 ']'
+ '[' 4925 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4925
+ '[' 4925 -lt 100000 ']'
+ FICHIER=velocity_Z_it04925
+ '[' 4925 -lt 10000 ']'
+ FICHIER=velocity_Z_it004925
+ '[' 4925 -lt 1000 ']'
+ '[' 4925 -lt 100 ']'
+ echo velocity_Z_it004925
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004925 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4925+25
+ '[' 4950 -le 8500 ']'
+ '[' 4950 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4950
+ '[' 4950 -lt 100000 ']'
+ FICHIER=velocity_Z_it04950
+ '[' 4950 -lt 10000 ']'
+ FICHIER=velocity_Z_it004950
+ '[' 4950 -lt 1000 ']'
+ '[' 4950 -lt 100 ']'
+ echo velocity_Z_it004950
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004950 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4950+25
+ '[' 4975 -le 8500 ']'
+ '[' 4975 -lt 1000000 ']'
+ FICHIER=velocity_Z_it4975
+ '[' 4975 -lt 100000 ']'
+ FICHIER=velocity_Z_it04975
+ '[' 4975 -lt 10000 ']'
+ FICHIER=velocity_Z_it004975
+ '[' 4975 -lt 1000 ']'
+ '[' 4975 -lt 100 ']'
+ echo velocity_Z_it004975
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it004975 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=4975+25
+ '[' 5000 -le 8500 ']'
+ '[' 5000 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5000
+ '[' 5000 -lt 100000 ']'
+ FICHIER=velocity_Z_it05000
+ '[' 5000 -lt 10000 ']'
+ FICHIER=velocity_Z_it005000
+ '[' 5000 -lt 1000 ']'
+ '[' 5000 -lt 100 ']'
+ echo velocity_Z_it005000
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005000 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5000+25
+ '[' 5025 -le 8500 ']'
+ '[' 5025 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5025
+ '[' 5025 -lt 100000 ']'
+ FICHIER=velocity_Z_it05025
+ '[' 5025 -lt 10000 ']'
+ FICHIER=velocity_Z_it005025
+ '[' 5025 -lt 1000 ']'
+ '[' 5025 -lt 100 ']'
+ echo velocity_Z_it005025
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005025 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5025+25
+ '[' 5050 -le 8500 ']'
+ '[' 5050 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5050
+ '[' 5050 -lt 100000 ']'
+ FICHIER=velocity_Z_it05050
+ '[' 5050 -lt 10000 ']'
+ FICHIER=velocity_Z_it005050
+ '[' 5050 -lt 1000 ']'
+ '[' 5050 -lt 100 ']'
+ echo velocity_Z_it005050
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005050 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5050+25
+ '[' 5075 -le 8500 ']'
+ '[' 5075 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5075
+ '[' 5075 -lt 100000 ']'
+ FICHIER=velocity_Z_it05075
+ '[' 5075 -lt 10000 ']'
+ FICHIER=velocity_Z_it005075
+ '[' 5075 -lt 1000 ']'
+ '[' 5075 -lt 100 ']'
+ echo velocity_Z_it005075
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005075 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5075+25
+ '[' 5100 -le 8500 ']'
+ '[' 5100 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5100
+ '[' 5100 -lt 100000 ']'
+ FICHIER=velocity_Z_it05100
+ '[' 5100 -lt 10000 ']'
+ FICHIER=velocity_Z_it005100
+ '[' 5100 -lt 1000 ']'
+ '[' 5100 -lt 100 ']'
+ echo velocity_Z_it005100
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005100 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5100+25
+ '[' 5125 -le 8500 ']'
+ '[' 5125 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5125
+ '[' 5125 -lt 100000 ']'
+ FICHIER=velocity_Z_it05125
+ '[' 5125 -lt 10000 ']'
+ FICHIER=velocity_Z_it005125
+ '[' 5125 -lt 1000 ']'
+ '[' 5125 -lt 100 ']'
+ echo velocity_Z_it005125
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005125 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5125+25
+ '[' 5150 -le 8500 ']'
+ '[' 5150 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5150
+ '[' 5150 -lt 100000 ']'
+ FICHIER=velocity_Z_it05150
+ '[' 5150 -lt 10000 ']'
+ FICHIER=velocity_Z_it005150
+ '[' 5150 -lt 1000 ']'
+ '[' 5150 -lt 100 ']'
+ echo velocity_Z_it005150
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005150 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5150+25
+ '[' 5175 -le 8500 ']'
+ '[' 5175 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5175
+ '[' 5175 -lt 100000 ']'
+ FICHIER=velocity_Z_it05175
+ '[' 5175 -lt 10000 ']'
+ FICHIER=velocity_Z_it005175
+ '[' 5175 -lt 1000 ']'
+ '[' 5175 -lt 100 ']'
+ echo velocity_Z_it005175
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005175 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5175+25
+ '[' 5200 -le 8500 ']'
+ '[' 5200 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5200
+ '[' 5200 -lt 100000 ']'
+ FICHIER=velocity_Z_it05200
+ '[' 5200 -lt 10000 ']'
+ FICHIER=velocity_Z_it005200
+ '[' 5200 -lt 1000 ']'
+ '[' 5200 -lt 100 ']'
+ echo velocity_Z_it005200
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005200 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5200+25
+ '[' 5225 -le 8500 ']'
+ '[' 5225 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5225
+ '[' 5225 -lt 100000 ']'
+ FICHIER=velocity_Z_it05225
+ '[' 5225 -lt 10000 ']'
+ FICHIER=velocity_Z_it005225
+ '[' 5225 -lt 1000 ']'
+ '[' 5225 -lt 100 ']'
+ echo velocity_Z_it005225
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005225 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5225+25
+ '[' 5250 -le 8500 ']'
+ '[' 5250 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5250
+ '[' 5250 -lt 100000 ']'
+ FICHIER=velocity_Z_it05250
+ '[' 5250 -lt 10000 ']'
+ FICHIER=velocity_Z_it005250
+ '[' 5250 -lt 1000 ']'
+ '[' 5250 -lt 100 ']'
+ echo velocity_Z_it005250
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005250 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5250+25
+ '[' 5275 -le 8500 ']'
+ '[' 5275 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5275
+ '[' 5275 -lt 100000 ']'
+ FICHIER=velocity_Z_it05275
+ '[' 5275 -lt 10000 ']'
+ FICHIER=velocity_Z_it005275
+ '[' 5275 -lt 1000 ']'
+ '[' 5275 -lt 100 ']'
+ echo velocity_Z_it005275
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005275 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5275+25
+ '[' 5300 -le 8500 ']'
+ '[' 5300 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5300
+ '[' 5300 -lt 100000 ']'
+ FICHIER=velocity_Z_it05300
+ '[' 5300 -lt 10000 ']'
+ FICHIER=velocity_Z_it005300
+ '[' 5300 -lt 1000 ']'
+ '[' 5300 -lt 100 ']'
+ echo velocity_Z_it005300
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005300 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5300+25
+ '[' 5325 -le 8500 ']'
+ '[' 5325 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5325
+ '[' 5325 -lt 100000 ']'
+ FICHIER=velocity_Z_it05325
+ '[' 5325 -lt 10000 ']'
+ FICHIER=velocity_Z_it005325
+ '[' 5325 -lt 1000 ']'
+ '[' 5325 -lt 100 ']'
+ echo velocity_Z_it005325
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005325 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5325+25
+ '[' 5350 -le 8500 ']'
+ '[' 5350 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5350
+ '[' 5350 -lt 100000 ']'
+ FICHIER=velocity_Z_it05350
+ '[' 5350 -lt 10000 ']'
+ FICHIER=velocity_Z_it005350
+ '[' 5350 -lt 1000 ']'
+ '[' 5350 -lt 100 ']'
+ echo velocity_Z_it005350
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005350 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5350+25
+ '[' 5375 -le 8500 ']'
+ '[' 5375 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5375
+ '[' 5375 -lt 100000 ']'
+ FICHIER=velocity_Z_it05375
+ '[' 5375 -lt 10000 ']'
+ FICHIER=velocity_Z_it005375
+ '[' 5375 -lt 1000 ']'
+ '[' 5375 -lt 100 ']'
+ echo velocity_Z_it005375
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005375 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5375+25
+ '[' 5400 -le 8500 ']'
+ '[' 5400 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5400
+ '[' 5400 -lt 100000 ']'
+ FICHIER=velocity_Z_it05400
+ '[' 5400 -lt 10000 ']'
+ FICHIER=velocity_Z_it005400
+ '[' 5400 -lt 1000 ']'
+ '[' 5400 -lt 100 ']'
+ echo velocity_Z_it005400
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005400 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5400+25
+ '[' 5425 -le 8500 ']'
+ '[' 5425 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5425
+ '[' 5425 -lt 100000 ']'
+ FICHIER=velocity_Z_it05425
+ '[' 5425 -lt 10000 ']'
+ FICHIER=velocity_Z_it005425
+ '[' 5425 -lt 1000 ']'
+ '[' 5425 -lt 100 ']'
+ echo velocity_Z_it005425
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005425 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5425+25
+ '[' 5450 -le 8500 ']'
+ '[' 5450 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5450
+ '[' 5450 -lt 100000 ']'
+ FICHIER=velocity_Z_it05450
+ '[' 5450 -lt 10000 ']'
+ FICHIER=velocity_Z_it005450
+ '[' 5450 -lt 1000 ']'
+ '[' 5450 -lt 100 ']'
+ echo velocity_Z_it005450
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005450 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5450+25
+ '[' 5475 -le 8500 ']'
+ '[' 5475 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5475
+ '[' 5475 -lt 100000 ']'
+ FICHIER=velocity_Z_it05475
+ '[' 5475 -lt 10000 ']'
+ FICHIER=velocity_Z_it005475
+ '[' 5475 -lt 1000 ']'
+ '[' 5475 -lt 100 ']'
+ echo velocity_Z_it005475
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005475 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5475+25
+ '[' 5500 -le 8500 ']'
+ '[' 5500 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5500
+ '[' 5500 -lt 100000 ']'
+ FICHIER=velocity_Z_it05500
+ '[' 5500 -lt 10000 ']'
+ FICHIER=velocity_Z_it005500
+ '[' 5500 -lt 1000 ']'
+ '[' 5500 -lt 100 ']'
+ echo velocity_Z_it005500
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005500 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5500+25
+ '[' 5525 -le 8500 ']'
+ '[' 5525 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5525
+ '[' 5525 -lt 100000 ']'
+ FICHIER=velocity_Z_it05525
+ '[' 5525 -lt 10000 ']'
+ FICHIER=velocity_Z_it005525
+ '[' 5525 -lt 1000 ']'
+ '[' 5525 -lt 100 ']'
+ echo velocity_Z_it005525
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005525 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5525+25
+ '[' 5550 -le 8500 ']'
+ '[' 5550 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5550
+ '[' 5550 -lt 100000 ']'
+ FICHIER=velocity_Z_it05550
+ '[' 5550 -lt 10000 ']'
+ FICHIER=velocity_Z_it005550
+ '[' 5550 -lt 1000 ']'
+ '[' 5550 -lt 100 ']'
+ echo velocity_Z_it005550
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005550 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5550+25
+ '[' 5575 -le 8500 ']'
+ '[' 5575 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5575
+ '[' 5575 -lt 100000 ']'
+ FICHIER=velocity_Z_it05575
+ '[' 5575 -lt 10000 ']'
+ FICHIER=velocity_Z_it005575
+ '[' 5575 -lt 1000 ']'
+ '[' 5575 -lt 100 ']'
+ echo velocity_Z_it005575
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005575 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5575+25
+ '[' 5600 -le 8500 ']'
+ '[' 5600 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5600
+ '[' 5600 -lt 100000 ']'
+ FICHIER=velocity_Z_it05600
+ '[' 5600 -lt 10000 ']'
+ FICHIER=velocity_Z_it005600
+ '[' 5600 -lt 1000 ']'
+ '[' 5600 -lt 100 ']'
+ echo velocity_Z_it005600
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005600 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5600+25
+ '[' 5625 -le 8500 ']'
+ '[' 5625 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5625
+ '[' 5625 -lt 100000 ']'
+ FICHIER=velocity_Z_it05625
+ '[' 5625 -lt 10000 ']'
+ FICHIER=velocity_Z_it005625
+ '[' 5625 -lt 1000 ']'
+ '[' 5625 -lt 100 ']'
+ echo velocity_Z_it005625
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005625 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5625+25
+ '[' 5650 -le 8500 ']'
+ '[' 5650 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5650
+ '[' 5650 -lt 100000 ']'
+ FICHIER=velocity_Z_it05650
+ '[' 5650 -lt 10000 ']'
+ FICHIER=velocity_Z_it005650
+ '[' 5650 -lt 1000 ']'
+ '[' 5650 -lt 100 ']'
+ echo velocity_Z_it005650
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005650 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5650+25
+ '[' 5675 -le 8500 ']'
+ '[' 5675 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5675
+ '[' 5675 -lt 100000 ']'
+ FICHIER=velocity_Z_it05675
+ '[' 5675 -lt 10000 ']'
+ FICHIER=velocity_Z_it005675
+ '[' 5675 -lt 1000 ']'
+ '[' 5675 -lt 100 ']'
+ echo velocity_Z_it005675
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005675 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5675+25
+ '[' 5700 -le 8500 ']'
+ '[' 5700 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5700
+ '[' 5700 -lt 100000 ']'
+ FICHIER=velocity_Z_it05700
+ '[' 5700 -lt 10000 ']'
+ FICHIER=velocity_Z_it005700
+ '[' 5700 -lt 1000 ']'
+ '[' 5700 -lt 100 ']'
+ echo velocity_Z_it005700
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005700 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5700+25
+ '[' 5725 -le 8500 ']'
+ '[' 5725 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5725
+ '[' 5725 -lt 100000 ']'
+ FICHIER=velocity_Z_it05725
+ '[' 5725 -lt 10000 ']'
+ FICHIER=velocity_Z_it005725
+ '[' 5725 -lt 1000 ']'
+ '[' 5725 -lt 100 ']'
+ echo velocity_Z_it005725
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005725 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5725+25
+ '[' 5750 -le 8500 ']'
+ '[' 5750 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5750
+ '[' 5750 -lt 100000 ']'
+ FICHIER=velocity_Z_it05750
+ '[' 5750 -lt 10000 ']'
+ FICHIER=velocity_Z_it005750
+ '[' 5750 -lt 1000 ']'
+ '[' 5750 -lt 100 ']'
+ echo velocity_Z_it005750
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005750 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5750+25
+ '[' 5775 -le 8500 ']'
+ '[' 5775 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5775
+ '[' 5775 -lt 100000 ']'
+ FICHIER=velocity_Z_it05775
+ '[' 5775 -lt 10000 ']'
+ FICHIER=velocity_Z_it005775
+ '[' 5775 -lt 1000 ']'
+ '[' 5775 -lt 100 ']'
+ echo velocity_Z_it005775
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005775 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5775+25
+ '[' 5800 -le 8500 ']'
+ '[' 5800 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5800
+ '[' 5800 -lt 100000 ']'
+ FICHIER=velocity_Z_it05800
+ '[' 5800 -lt 10000 ']'
+ FICHIER=velocity_Z_it005800
+ '[' 5800 -lt 1000 ']'
+ '[' 5800 -lt 100 ']'
+ echo velocity_Z_it005800
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005800 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5800+25
+ '[' 5825 -le 8500 ']'
+ '[' 5825 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5825
+ '[' 5825 -lt 100000 ']'
+ FICHIER=velocity_Z_it05825
+ '[' 5825 -lt 10000 ']'
+ FICHIER=velocity_Z_it005825
+ '[' 5825 -lt 1000 ']'
+ '[' 5825 -lt 100 ']'
+ echo velocity_Z_it005825
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005825 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5825+25
+ '[' 5850 -le 8500 ']'
+ '[' 5850 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5850
+ '[' 5850 -lt 100000 ']'
+ FICHIER=velocity_Z_it05850
+ '[' 5850 -lt 10000 ']'
+ FICHIER=velocity_Z_it005850
+ '[' 5850 -lt 1000 ']'
+ '[' 5850 -lt 100 ']'
+ echo velocity_Z_it005850
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005850 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5850+25
+ '[' 5875 -le 8500 ']'
+ '[' 5875 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5875
+ '[' 5875 -lt 100000 ']'
+ FICHIER=velocity_Z_it05875
+ '[' 5875 -lt 10000 ']'
+ FICHIER=velocity_Z_it005875
+ '[' 5875 -lt 1000 ']'
+ '[' 5875 -lt 100 ']'
+ echo velocity_Z_it005875
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005875 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5875+25
+ '[' 5900 -le 8500 ']'
+ '[' 5900 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5900
+ '[' 5900 -lt 100000 ']'
+ FICHIER=velocity_Z_it05900
+ '[' 5900 -lt 10000 ']'
+ FICHIER=velocity_Z_it005900
+ '[' 5900 -lt 1000 ']'
+ '[' 5900 -lt 100 ']'
+ echo velocity_Z_it005900
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005900 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5900+25
+ '[' 5925 -le 8500 ']'
+ '[' 5925 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5925
+ '[' 5925 -lt 100000 ']'
+ FICHIER=velocity_Z_it05925
+ '[' 5925 -lt 10000 ']'
+ FICHIER=velocity_Z_it005925
+ '[' 5925 -lt 1000 ']'
+ '[' 5925 -lt 100 ']'
+ echo velocity_Z_it005925
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005925 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5925+25
+ '[' 5950 -le 8500 ']'
+ '[' 5950 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5950
+ '[' 5950 -lt 100000 ']'
+ FICHIER=velocity_Z_it05950
+ '[' 5950 -lt 10000 ']'
+ FICHIER=velocity_Z_it005950
+ '[' 5950 -lt 1000 ']'
+ '[' 5950 -lt 100 ']'
+ echo velocity_Z_it005950
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005950 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5950+25
+ '[' 5975 -le 8500 ']'
+ '[' 5975 -lt 1000000 ']'
+ FICHIER=velocity_Z_it5975
+ '[' 5975 -lt 100000 ']'
+ FICHIER=velocity_Z_it05975
+ '[' 5975 -lt 10000 ']'
+ FICHIER=velocity_Z_it005975
+ '[' 5975 -lt 1000 ']'
+ '[' 5975 -lt 100 ']'
+ echo velocity_Z_it005975
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it005975 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=5975+25
+ '[' 6000 -le 8500 ']'
+ '[' 6000 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6000
+ '[' 6000 -lt 100000 ']'
+ FICHIER=velocity_Z_it06000
+ '[' 6000 -lt 10000 ']'
+ FICHIER=velocity_Z_it006000
+ '[' 6000 -lt 1000 ']'
+ '[' 6000 -lt 100 ']'
+ echo velocity_Z_it006000
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006000 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6000+25
+ '[' 6025 -le 8500 ']'
+ '[' 6025 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6025
+ '[' 6025 -lt 100000 ']'
+ FICHIER=velocity_Z_it06025
+ '[' 6025 -lt 10000 ']'
+ FICHIER=velocity_Z_it006025
+ '[' 6025 -lt 1000 ']'
+ '[' 6025 -lt 100 ']'
+ echo velocity_Z_it006025
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006025 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6025+25
+ '[' 6050 -le 8500 ']'
+ '[' 6050 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6050
+ '[' 6050 -lt 100000 ']'
+ FICHIER=velocity_Z_it06050
+ '[' 6050 -lt 10000 ']'
+ FICHIER=velocity_Z_it006050
+ '[' 6050 -lt 1000 ']'
+ '[' 6050 -lt 100 ']'
+ echo velocity_Z_it006050
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006050 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6050+25
+ '[' 6075 -le 8500 ']'
+ '[' 6075 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6075
+ '[' 6075 -lt 100000 ']'
+ FICHIER=velocity_Z_it06075
+ '[' 6075 -lt 10000 ']'
+ FICHIER=velocity_Z_it006075
+ '[' 6075 -lt 1000 ']'
+ '[' 6075 -lt 100 ']'
+ echo velocity_Z_it006075
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006075 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6075+25
+ '[' 6100 -le 8500 ']'
+ '[' 6100 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6100
+ '[' 6100 -lt 100000 ']'
+ FICHIER=velocity_Z_it06100
+ '[' 6100 -lt 10000 ']'
+ FICHIER=velocity_Z_it006100
+ '[' 6100 -lt 1000 ']'
+ '[' 6100 -lt 100 ']'
+ echo velocity_Z_it006100
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006100 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6100+25
+ '[' 6125 -le 8500 ']'
+ '[' 6125 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6125
+ '[' 6125 -lt 100000 ']'
+ FICHIER=velocity_Z_it06125
+ '[' 6125 -lt 10000 ']'
+ FICHIER=velocity_Z_it006125
+ '[' 6125 -lt 1000 ']'
+ '[' 6125 -lt 100 ']'
+ echo velocity_Z_it006125
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006125 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6125+25
+ '[' 6150 -le 8500 ']'
+ '[' 6150 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6150
+ '[' 6150 -lt 100000 ']'
+ FICHIER=velocity_Z_it06150
+ '[' 6150 -lt 10000 ']'
+ FICHIER=velocity_Z_it006150
+ '[' 6150 -lt 1000 ']'
+ '[' 6150 -lt 100 ']'
+ echo velocity_Z_it006150
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006150 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6150+25
+ '[' 6175 -le 8500 ']'
+ '[' 6175 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6175
+ '[' 6175 -lt 100000 ']'
+ FICHIER=velocity_Z_it06175
+ '[' 6175 -lt 10000 ']'
+ FICHIER=velocity_Z_it006175
+ '[' 6175 -lt 1000 ']'
+ '[' 6175 -lt 100 ']'
+ echo velocity_Z_it006175
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006175 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6175+25
+ '[' 6200 -le 8500 ']'
+ '[' 6200 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6200
+ '[' 6200 -lt 100000 ']'
+ FICHIER=velocity_Z_it06200
+ '[' 6200 -lt 10000 ']'
+ FICHIER=velocity_Z_it006200
+ '[' 6200 -lt 1000 ']'
+ '[' 6200 -lt 100 ']'
+ echo velocity_Z_it006200
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006200 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6200+25
+ '[' 6225 -le 8500 ']'
+ '[' 6225 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6225
+ '[' 6225 -lt 100000 ']'
+ FICHIER=velocity_Z_it06225
+ '[' 6225 -lt 10000 ']'
+ FICHIER=velocity_Z_it006225
+ '[' 6225 -lt 1000 ']'
+ '[' 6225 -lt 100 ']'
+ echo velocity_Z_it006225
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006225 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6225+25
+ '[' 6250 -le 8500 ']'
+ '[' 6250 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6250
+ '[' 6250 -lt 100000 ']'
+ FICHIER=velocity_Z_it06250
+ '[' 6250 -lt 10000 ']'
+ FICHIER=velocity_Z_it006250
+ '[' 6250 -lt 1000 ']'
+ '[' 6250 -lt 100 ']'
+ echo velocity_Z_it006250
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006250 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6250+25
+ '[' 6275 -le 8500 ']'
+ '[' 6275 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6275
+ '[' 6275 -lt 100000 ']'
+ FICHIER=velocity_Z_it06275
+ '[' 6275 -lt 10000 ']'
+ FICHIER=velocity_Z_it006275
+ '[' 6275 -lt 1000 ']'
+ '[' 6275 -lt 100 ']'
+ echo velocity_Z_it006275
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006275 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6275+25
+ '[' 6300 -le 8500 ']'
+ '[' 6300 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6300
+ '[' 6300 -lt 100000 ']'
+ FICHIER=velocity_Z_it06300
+ '[' 6300 -lt 10000 ']'
+ FICHIER=velocity_Z_it006300
+ '[' 6300 -lt 1000 ']'
+ '[' 6300 -lt 100 ']'
+ echo velocity_Z_it006300
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006300 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6300+25
+ '[' 6325 -le 8500 ']'
+ '[' 6325 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6325
+ '[' 6325 -lt 100000 ']'
+ FICHIER=velocity_Z_it06325
+ '[' 6325 -lt 10000 ']'
+ FICHIER=velocity_Z_it006325
+ '[' 6325 -lt 1000 ']'
+ '[' 6325 -lt 100 ']'
+ echo velocity_Z_it006325
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006325 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6325+25
+ '[' 6350 -le 8500 ']'
+ '[' 6350 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6350
+ '[' 6350 -lt 100000 ']'
+ FICHIER=velocity_Z_it06350
+ '[' 6350 -lt 10000 ']'
+ FICHIER=velocity_Z_it006350
+ '[' 6350 -lt 1000 ']'
+ '[' 6350 -lt 100 ']'
+ echo velocity_Z_it006350
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006350 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6350+25
+ '[' 6375 -le 8500 ']'
+ '[' 6375 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6375
+ '[' 6375 -lt 100000 ']'
+ FICHIER=velocity_Z_it06375
+ '[' 6375 -lt 10000 ']'
+ FICHIER=velocity_Z_it006375
+ '[' 6375 -lt 1000 ']'
+ '[' 6375 -lt 100 ']'
+ echo velocity_Z_it006375
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006375 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6375+25
+ '[' 6400 -le 8500 ']'
+ '[' 6400 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6400
+ '[' 6400 -lt 100000 ']'
+ FICHIER=velocity_Z_it06400
+ '[' 6400 -lt 10000 ']'
+ FICHIER=velocity_Z_it006400
+ '[' 6400 -lt 1000 ']'
+ '[' 6400 -lt 100 ']'
+ echo velocity_Z_it006400
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006400 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6400+25
+ '[' 6425 -le 8500 ']'
+ '[' 6425 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6425
+ '[' 6425 -lt 100000 ']'
+ FICHIER=velocity_Z_it06425
+ '[' 6425 -lt 10000 ']'
+ FICHIER=velocity_Z_it006425
+ '[' 6425 -lt 1000 ']'
+ '[' 6425 -lt 100 ']'
+ echo velocity_Z_it006425
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006425 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6425+25
+ '[' 6450 -le 8500 ']'
+ '[' 6450 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6450
+ '[' 6450 -lt 100000 ']'
+ FICHIER=velocity_Z_it06450
+ '[' 6450 -lt 10000 ']'
+ FICHIER=velocity_Z_it006450
+ '[' 6450 -lt 1000 ']'
+ '[' 6450 -lt 100 ']'
+ echo velocity_Z_it006450
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006450 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6450+25
+ '[' 6475 -le 8500 ']'
+ '[' 6475 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6475
+ '[' 6475 -lt 100000 ']'
+ FICHIER=velocity_Z_it06475
+ '[' 6475 -lt 10000 ']'
+ FICHIER=velocity_Z_it006475
+ '[' 6475 -lt 1000 ']'
+ '[' 6475 -lt 100 ']'
+ echo velocity_Z_it006475
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006475 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6475+25
+ '[' 6500 -le 8500 ']'
+ '[' 6500 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6500
+ '[' 6500 -lt 100000 ']'
+ FICHIER=velocity_Z_it06500
+ '[' 6500 -lt 10000 ']'
+ FICHIER=velocity_Z_it006500
+ '[' 6500 -lt 1000 ']'
+ '[' 6500 -lt 100 ']'
+ echo velocity_Z_it006500
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006500 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6500+25
+ '[' 6525 -le 8500 ']'
+ '[' 6525 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6525
+ '[' 6525 -lt 100000 ']'
+ FICHIER=velocity_Z_it06525
+ '[' 6525 -lt 10000 ']'
+ FICHIER=velocity_Z_it006525
+ '[' 6525 -lt 1000 ']'
+ '[' 6525 -lt 100 ']'
+ echo velocity_Z_it006525
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006525 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6525+25
+ '[' 6550 -le 8500 ']'
+ '[' 6550 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6550
+ '[' 6550 -lt 100000 ']'
+ FICHIER=velocity_Z_it06550
+ '[' 6550 -lt 10000 ']'
+ FICHIER=velocity_Z_it006550
+ '[' 6550 -lt 1000 ']'
+ '[' 6550 -lt 100 ']'
+ echo velocity_Z_it006550
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006550 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6550+25
+ '[' 6575 -le 8500 ']'
+ '[' 6575 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6575
+ '[' 6575 -lt 100000 ']'
+ FICHIER=velocity_Z_it06575
+ '[' 6575 -lt 10000 ']'
+ FICHIER=velocity_Z_it006575
+ '[' 6575 -lt 1000 ']'
+ '[' 6575 -lt 100 ']'
+ echo velocity_Z_it006575
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006575 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6575+25
+ '[' 6600 -le 8500 ']'
+ '[' 6600 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6600
+ '[' 6600 -lt 100000 ']'
+ FICHIER=velocity_Z_it06600
+ '[' 6600 -lt 10000 ']'
+ FICHIER=velocity_Z_it006600
+ '[' 6600 -lt 1000 ']'
+ '[' 6600 -lt 100 ']'
+ echo velocity_Z_it006600
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006600 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6600+25
+ '[' 6625 -le 8500 ']'
+ '[' 6625 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6625
+ '[' 6625 -lt 100000 ']'
+ FICHIER=velocity_Z_it06625
+ '[' 6625 -lt 10000 ']'
+ FICHIER=velocity_Z_it006625
+ '[' 6625 -lt 1000 ']'
+ '[' 6625 -lt 100 ']'
+ echo velocity_Z_it006625
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006625 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6625+25
+ '[' 6650 -le 8500 ']'
+ '[' 6650 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6650
+ '[' 6650 -lt 100000 ']'
+ FICHIER=velocity_Z_it06650
+ '[' 6650 -lt 10000 ']'
+ FICHIER=velocity_Z_it006650
+ '[' 6650 -lt 1000 ']'
+ '[' 6650 -lt 100 ']'
+ echo velocity_Z_it006650
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006650 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6650+25
+ '[' 6675 -le 8500 ']'
+ '[' 6675 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6675
+ '[' 6675 -lt 100000 ']'
+ FICHIER=velocity_Z_it06675
+ '[' 6675 -lt 10000 ']'
+ FICHIER=velocity_Z_it006675
+ '[' 6675 -lt 1000 ']'
+ '[' 6675 -lt 100 ']'
+ echo velocity_Z_it006675
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006675 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6675+25
+ '[' 6700 -le 8500 ']'
+ '[' 6700 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6700
+ '[' 6700 -lt 100000 ']'
+ FICHIER=velocity_Z_it06700
+ '[' 6700 -lt 10000 ']'
+ FICHIER=velocity_Z_it006700
+ '[' 6700 -lt 1000 ']'
+ '[' 6700 -lt 100 ']'
+ echo velocity_Z_it006700
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006700 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6700+25
+ '[' 6725 -le 8500 ']'
+ '[' 6725 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6725
+ '[' 6725 -lt 100000 ']'
+ FICHIER=velocity_Z_it06725
+ '[' 6725 -lt 10000 ']'
+ FICHIER=velocity_Z_it006725
+ '[' 6725 -lt 1000 ']'
+ '[' 6725 -lt 100 ']'
+ echo velocity_Z_it006725
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006725 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6725+25
+ '[' 6750 -le 8500 ']'
+ '[' 6750 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6750
+ '[' 6750 -lt 100000 ']'
+ FICHIER=velocity_Z_it06750
+ '[' 6750 -lt 10000 ']'
+ FICHIER=velocity_Z_it006750
+ '[' 6750 -lt 1000 ']'
+ '[' 6750 -lt 100 ']'
+ echo velocity_Z_it006750
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006750 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6750+25
+ '[' 6775 -le 8500 ']'
+ '[' 6775 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6775
+ '[' 6775 -lt 100000 ']'
+ FICHIER=velocity_Z_it06775
+ '[' 6775 -lt 10000 ']'
+ FICHIER=velocity_Z_it006775
+ '[' 6775 -lt 1000 ']'
+ '[' 6775 -lt 100 ']'
+ echo velocity_Z_it006775
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006775 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6775+25
+ '[' 6800 -le 8500 ']'
+ '[' 6800 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6800
+ '[' 6800 -lt 100000 ']'
+ FICHIER=velocity_Z_it06800
+ '[' 6800 -lt 10000 ']'
+ FICHIER=velocity_Z_it006800
+ '[' 6800 -lt 1000 ']'
+ '[' 6800 -lt 100 ']'
+ echo velocity_Z_it006800
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006800 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6800+25
+ '[' 6825 -le 8500 ']'
+ '[' 6825 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6825
+ '[' 6825 -lt 100000 ']'
+ FICHIER=velocity_Z_it06825
+ '[' 6825 -lt 10000 ']'
+ FICHIER=velocity_Z_it006825
+ '[' 6825 -lt 1000 ']'
+ '[' 6825 -lt 100 ']'
+ echo velocity_Z_it006825
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006825 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6825+25
+ '[' 6850 -le 8500 ']'
+ '[' 6850 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6850
+ '[' 6850 -lt 100000 ']'
+ FICHIER=velocity_Z_it06850
+ '[' 6850 -lt 10000 ']'
+ FICHIER=velocity_Z_it006850
+ '[' 6850 -lt 1000 ']'
+ '[' 6850 -lt 100 ']'
+ echo velocity_Z_it006850
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006850 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6850+25
+ '[' 6875 -le 8500 ']'
+ '[' 6875 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6875
+ '[' 6875 -lt 100000 ']'
+ FICHIER=velocity_Z_it06875
+ '[' 6875 -lt 10000 ']'
+ FICHIER=velocity_Z_it006875
+ '[' 6875 -lt 1000 ']'
+ '[' 6875 -lt 100 ']'
+ echo velocity_Z_it006875
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006875 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6875+25
+ '[' 6900 -le 8500 ']'
+ '[' 6900 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6900
+ '[' 6900 -lt 100000 ']'
+ FICHIER=velocity_Z_it06900
+ '[' 6900 -lt 10000 ']'
+ FICHIER=velocity_Z_it006900
+ '[' 6900 -lt 1000 ']'
+ '[' 6900 -lt 100 ']'
+ echo velocity_Z_it006900
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006900 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6900+25
+ '[' 6925 -le 8500 ']'
+ '[' 6925 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6925
+ '[' 6925 -lt 100000 ']'
+ FICHIER=velocity_Z_it06925
+ '[' 6925 -lt 10000 ']'
+ FICHIER=velocity_Z_it006925
+ '[' 6925 -lt 1000 ']'
+ '[' 6925 -lt 100 ']'
+ echo velocity_Z_it006925
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006925 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6925+25
+ '[' 6950 -le 8500 ']'
+ '[' 6950 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6950
+ '[' 6950 -lt 100000 ']'
+ FICHIER=velocity_Z_it06950
+ '[' 6950 -lt 10000 ']'
+ FICHIER=velocity_Z_it006950
+ '[' 6950 -lt 1000 ']'
+ '[' 6950 -lt 100 ']'
+ echo velocity_Z_it006950
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006950 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6950+25
+ '[' 6975 -le 8500 ']'
+ '[' 6975 -lt 1000000 ']'
+ FICHIER=velocity_Z_it6975
+ '[' 6975 -lt 100000 ']'
+ FICHIER=velocity_Z_it06975
+ '[' 6975 -lt 10000 ']'
+ FICHIER=velocity_Z_it006975
+ '[' 6975 -lt 1000 ']'
+ '[' 6975 -lt 100 ']'
+ echo velocity_Z_it006975
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it006975 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=6975+25
+ '[' 7000 -le 8500 ']'
+ '[' 7000 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7000
+ '[' 7000 -lt 100000 ']'
+ FICHIER=velocity_Z_it07000
+ '[' 7000 -lt 10000 ']'
+ FICHIER=velocity_Z_it007000
+ '[' 7000 -lt 1000 ']'
+ '[' 7000 -lt 100 ']'
+ echo velocity_Z_it007000
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007000 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7000+25
+ '[' 7025 -le 8500 ']'
+ '[' 7025 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7025
+ '[' 7025 -lt 100000 ']'
+ FICHIER=velocity_Z_it07025
+ '[' 7025 -lt 10000 ']'
+ FICHIER=velocity_Z_it007025
+ '[' 7025 -lt 1000 ']'
+ '[' 7025 -lt 100 ']'
+ echo velocity_Z_it007025
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007025 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7025+25
+ '[' 7050 -le 8500 ']'
+ '[' 7050 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7050
+ '[' 7050 -lt 100000 ']'
+ FICHIER=velocity_Z_it07050
+ '[' 7050 -lt 10000 ']'
+ FICHIER=velocity_Z_it007050
+ '[' 7050 -lt 1000 ']'
+ '[' 7050 -lt 100 ']'
+ echo velocity_Z_it007050
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007050 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7050+25
+ '[' 7075 -le 8500 ']'
+ '[' 7075 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7075
+ '[' 7075 -lt 100000 ']'
+ FICHIER=velocity_Z_it07075
+ '[' 7075 -lt 10000 ']'
+ FICHIER=velocity_Z_it007075
+ '[' 7075 -lt 1000 ']'
+ '[' 7075 -lt 100 ']'
+ echo velocity_Z_it007075
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007075 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7075+25
+ '[' 7100 -le 8500 ']'
+ '[' 7100 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7100
+ '[' 7100 -lt 100000 ']'
+ FICHIER=velocity_Z_it07100
+ '[' 7100 -lt 10000 ']'
+ FICHIER=velocity_Z_it007100
+ '[' 7100 -lt 1000 ']'
+ '[' 7100 -lt 100 ']'
+ echo velocity_Z_it007100
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007100 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7100+25
+ '[' 7125 -le 8500 ']'
+ '[' 7125 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7125
+ '[' 7125 -lt 100000 ']'
+ FICHIER=velocity_Z_it07125
+ '[' 7125 -lt 10000 ']'
+ FICHIER=velocity_Z_it007125
+ '[' 7125 -lt 1000 ']'
+ '[' 7125 -lt 100 ']'
+ echo velocity_Z_it007125
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007125 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7125+25
+ '[' 7150 -le 8500 ']'
+ '[' 7150 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7150
+ '[' 7150 -lt 100000 ']'
+ FICHIER=velocity_Z_it07150
+ '[' 7150 -lt 10000 ']'
+ FICHIER=velocity_Z_it007150
+ '[' 7150 -lt 1000 ']'
+ '[' 7150 -lt 100 ']'
+ echo velocity_Z_it007150
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007150 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7150+25
+ '[' 7175 -le 8500 ']'
+ '[' 7175 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7175
+ '[' 7175 -lt 100000 ']'
+ FICHIER=velocity_Z_it07175
+ '[' 7175 -lt 10000 ']'
+ FICHIER=velocity_Z_it007175
+ '[' 7175 -lt 1000 ']'
+ '[' 7175 -lt 100 ']'
+ echo velocity_Z_it007175
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007175 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7175+25
+ '[' 7200 -le 8500 ']'
+ '[' 7200 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7200
+ '[' 7200 -lt 100000 ']'
+ FICHIER=velocity_Z_it07200
+ '[' 7200 -lt 10000 ']'
+ FICHIER=velocity_Z_it007200
+ '[' 7200 -lt 1000 ']'
+ '[' 7200 -lt 100 ']'
+ echo velocity_Z_it007200
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007200 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7200+25
+ '[' 7225 -le 8500 ']'
+ '[' 7225 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7225
+ '[' 7225 -lt 100000 ']'
+ FICHIER=velocity_Z_it07225
+ '[' 7225 -lt 10000 ']'
+ FICHIER=velocity_Z_it007225
+ '[' 7225 -lt 1000 ']'
+ '[' 7225 -lt 100 ']'
+ echo velocity_Z_it007225
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007225 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7225+25
+ '[' 7250 -le 8500 ']'
+ '[' 7250 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7250
+ '[' 7250 -lt 100000 ']'
+ FICHIER=velocity_Z_it07250
+ '[' 7250 -lt 10000 ']'
+ FICHIER=velocity_Z_it007250
+ '[' 7250 -lt 1000 ']'
+ '[' 7250 -lt 100 ']'
+ echo velocity_Z_it007250
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007250 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7250+25
+ '[' 7275 -le 8500 ']'
+ '[' 7275 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7275
+ '[' 7275 -lt 100000 ']'
+ FICHIER=velocity_Z_it07275
+ '[' 7275 -lt 10000 ']'
+ FICHIER=velocity_Z_it007275
+ '[' 7275 -lt 1000 ']'
+ '[' 7275 -lt 100 ']'
+ echo velocity_Z_it007275
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007275 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7275+25
+ '[' 7300 -le 8500 ']'
+ '[' 7300 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7300
+ '[' 7300 -lt 100000 ']'
+ FICHIER=velocity_Z_it07300
+ '[' 7300 -lt 10000 ']'
+ FICHIER=velocity_Z_it007300
+ '[' 7300 -lt 1000 ']'
+ '[' 7300 -lt 100 ']'
+ echo velocity_Z_it007300
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007300 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7300+25
+ '[' 7325 -le 8500 ']'
+ '[' 7325 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7325
+ '[' 7325 -lt 100000 ']'
+ FICHIER=velocity_Z_it07325
+ '[' 7325 -lt 10000 ']'
+ FICHIER=velocity_Z_it007325
+ '[' 7325 -lt 1000 ']'
+ '[' 7325 -lt 100 ']'
+ echo velocity_Z_it007325
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007325 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7325+25
+ '[' 7350 -le 8500 ']'
+ '[' 7350 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7350
+ '[' 7350 -lt 100000 ']'
+ FICHIER=velocity_Z_it07350
+ '[' 7350 -lt 10000 ']'
+ FICHIER=velocity_Z_it007350
+ '[' 7350 -lt 1000 ']'
+ '[' 7350 -lt 100 ']'
+ echo velocity_Z_it007350
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007350 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7350+25
+ '[' 7375 -le 8500 ']'
+ '[' 7375 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7375
+ '[' 7375 -lt 100000 ']'
+ FICHIER=velocity_Z_it07375
+ '[' 7375 -lt 10000 ']'
+ FICHIER=velocity_Z_it007375
+ '[' 7375 -lt 1000 ']'
+ '[' 7375 -lt 100 ']'
+ echo velocity_Z_it007375
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007375 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7375+25
+ '[' 7400 -le 8500 ']'
+ '[' 7400 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7400
+ '[' 7400 -lt 100000 ']'
+ FICHIER=velocity_Z_it07400
+ '[' 7400 -lt 10000 ']'
+ FICHIER=velocity_Z_it007400
+ '[' 7400 -lt 1000 ']'
+ '[' 7400 -lt 100 ']'
+ echo velocity_Z_it007400
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007400 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7400+25
+ '[' 7425 -le 8500 ']'
+ '[' 7425 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7425
+ '[' 7425 -lt 100000 ']'
+ FICHIER=velocity_Z_it07425
+ '[' 7425 -lt 10000 ']'
+ FICHIER=velocity_Z_it007425
+ '[' 7425 -lt 1000 ']'
+ '[' 7425 -lt 100 ']'
+ echo velocity_Z_it007425
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007425 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7425+25
+ '[' 7450 -le 8500 ']'
+ '[' 7450 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7450
+ '[' 7450 -lt 100000 ']'
+ FICHIER=velocity_Z_it07450
+ '[' 7450 -lt 10000 ']'
+ FICHIER=velocity_Z_it007450
+ '[' 7450 -lt 1000 ']'
+ '[' 7450 -lt 100 ']'
+ echo velocity_Z_it007450
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007450 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7450+25
+ '[' 7475 -le 8500 ']'
+ '[' 7475 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7475
+ '[' 7475 -lt 100000 ']'
+ FICHIER=velocity_Z_it07475
+ '[' 7475 -lt 10000 ']'
+ FICHIER=velocity_Z_it007475
+ '[' 7475 -lt 1000 ']'
+ '[' 7475 -lt 100 ']'
+ echo velocity_Z_it007475
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007475 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7475+25
+ '[' 7500 -le 8500 ']'
+ '[' 7500 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7500
+ '[' 7500 -lt 100000 ']'
+ FICHIER=velocity_Z_it07500
+ '[' 7500 -lt 10000 ']'
+ FICHIER=velocity_Z_it007500
+ '[' 7500 -lt 1000 ']'
+ '[' 7500 -lt 100 ']'
+ echo velocity_Z_it007500
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007500 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7500+25
+ '[' 7525 -le 8500 ']'
+ '[' 7525 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7525
+ '[' 7525 -lt 100000 ']'
+ FICHIER=velocity_Z_it07525
+ '[' 7525 -lt 10000 ']'
+ FICHIER=velocity_Z_it007525
+ '[' 7525 -lt 1000 ']'
+ '[' 7525 -lt 100 ']'
+ echo velocity_Z_it007525
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007525 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7525+25
+ '[' 7550 -le 8500 ']'
+ '[' 7550 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7550
+ '[' 7550 -lt 100000 ']'
+ FICHIER=velocity_Z_it07550
+ '[' 7550 -lt 10000 ']'
+ FICHIER=velocity_Z_it007550
+ '[' 7550 -lt 1000 ']'
+ '[' 7550 -lt 100 ']'
+ echo velocity_Z_it007550
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007550 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7550+25
+ '[' 7575 -le 8500 ']'
+ '[' 7575 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7575
+ '[' 7575 -lt 100000 ']'
+ FICHIER=velocity_Z_it07575
+ '[' 7575 -lt 10000 ']'
+ FICHIER=velocity_Z_it007575
+ '[' 7575 -lt 1000 ']'
+ '[' 7575 -lt 100 ']'
+ echo velocity_Z_it007575
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007575 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7575+25
+ '[' 7600 -le 8500 ']'
+ '[' 7600 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7600
+ '[' 7600 -lt 100000 ']'
+ FICHIER=velocity_Z_it07600
+ '[' 7600 -lt 10000 ']'
+ FICHIER=velocity_Z_it007600
+ '[' 7600 -lt 1000 ']'
+ '[' 7600 -lt 100 ']'
+ echo velocity_Z_it007600
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007600 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7600+25
+ '[' 7625 -le 8500 ']'
+ '[' 7625 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7625
+ '[' 7625 -lt 100000 ']'
+ FICHIER=velocity_Z_it07625
+ '[' 7625 -lt 10000 ']'
+ FICHIER=velocity_Z_it007625
+ '[' 7625 -lt 1000 ']'
+ '[' 7625 -lt 100 ']'
+ echo velocity_Z_it007625
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007625 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7625+25
+ '[' 7650 -le 8500 ']'
+ '[' 7650 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7650
+ '[' 7650 -lt 100000 ']'
+ FICHIER=velocity_Z_it07650
+ '[' 7650 -lt 10000 ']'
+ FICHIER=velocity_Z_it007650
+ '[' 7650 -lt 1000 ']'
+ '[' 7650 -lt 100 ']'
+ echo velocity_Z_it007650
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007650 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7650+25
+ '[' 7675 -le 8500 ']'
+ '[' 7675 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7675
+ '[' 7675 -lt 100000 ']'
+ FICHIER=velocity_Z_it07675
+ '[' 7675 -lt 10000 ']'
+ FICHIER=velocity_Z_it007675
+ '[' 7675 -lt 1000 ']'
+ '[' 7675 -lt 100 ']'
+ echo velocity_Z_it007675
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007675 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7675+25
+ '[' 7700 -le 8500 ']'
+ '[' 7700 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7700
+ '[' 7700 -lt 100000 ']'
+ FICHIER=velocity_Z_it07700
+ '[' 7700 -lt 10000 ']'
+ FICHIER=velocity_Z_it007700
+ '[' 7700 -lt 1000 ']'
+ '[' 7700 -lt 100 ']'
+ echo velocity_Z_it007700
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007700 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7700+25
+ '[' 7725 -le 8500 ']'
+ '[' 7725 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7725
+ '[' 7725 -lt 100000 ']'
+ FICHIER=velocity_Z_it07725
+ '[' 7725 -lt 10000 ']'
+ FICHIER=velocity_Z_it007725
+ '[' 7725 -lt 1000 ']'
+ '[' 7725 -lt 100 ']'
+ echo velocity_Z_it007725
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007725 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7725+25
+ '[' 7750 -le 8500 ']'
+ '[' 7750 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7750
+ '[' 7750 -lt 100000 ']'
+ FICHIER=velocity_Z_it07750
+ '[' 7750 -lt 10000 ']'
+ FICHIER=velocity_Z_it007750
+ '[' 7750 -lt 1000 ']'
+ '[' 7750 -lt 100 ']'
+ echo velocity_Z_it007750
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007750 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7750+25
+ '[' 7775 -le 8500 ']'
+ '[' 7775 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7775
+ '[' 7775 -lt 100000 ']'
+ FICHIER=velocity_Z_it07775
+ '[' 7775 -lt 10000 ']'
+ FICHIER=velocity_Z_it007775
+ '[' 7775 -lt 1000 ']'
+ '[' 7775 -lt 100 ']'
+ echo velocity_Z_it007775
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007775 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7775+25
+ '[' 7800 -le 8500 ']'
+ '[' 7800 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7800
+ '[' 7800 -lt 100000 ']'
+ FICHIER=velocity_Z_it07800
+ '[' 7800 -lt 10000 ']'
+ FICHIER=velocity_Z_it007800
+ '[' 7800 -lt 1000 ']'
+ '[' 7800 -lt 100 ']'
+ echo velocity_Z_it007800
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007800 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7800+25
+ '[' 7825 -le 8500 ']'
+ '[' 7825 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7825
+ '[' 7825 -lt 100000 ']'
+ FICHIER=velocity_Z_it07825
+ '[' 7825 -lt 10000 ']'
+ FICHIER=velocity_Z_it007825
+ '[' 7825 -lt 1000 ']'
+ '[' 7825 -lt 100 ']'
+ echo velocity_Z_it007825
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007825 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7825+25
+ '[' 7850 -le 8500 ']'
+ '[' 7850 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7850
+ '[' 7850 -lt 100000 ']'
+ FICHIER=velocity_Z_it07850
+ '[' 7850 -lt 10000 ']'
+ FICHIER=velocity_Z_it007850
+ '[' 7850 -lt 1000 ']'
+ '[' 7850 -lt 100 ']'
+ echo velocity_Z_it007850
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007850 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7850+25
+ '[' 7875 -le 8500 ']'
+ '[' 7875 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7875
+ '[' 7875 -lt 100000 ']'
+ FICHIER=velocity_Z_it07875
+ '[' 7875 -lt 10000 ']'
+ FICHIER=velocity_Z_it007875
+ '[' 7875 -lt 1000 ']'
+ '[' 7875 -lt 100 ']'
+ echo velocity_Z_it007875
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007875 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7875+25
+ '[' 7900 -le 8500 ']'
+ '[' 7900 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7900
+ '[' 7900 -lt 100000 ']'
+ FICHIER=velocity_Z_it07900
+ '[' 7900 -lt 10000 ']'
+ FICHIER=velocity_Z_it007900
+ '[' 7900 -lt 1000 ']'
+ '[' 7900 -lt 100 ']'
+ echo velocity_Z_it007900
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007900 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7900+25
+ '[' 7925 -le 8500 ']'
+ '[' 7925 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7925
+ '[' 7925 -lt 100000 ']'
+ FICHIER=velocity_Z_it07925
+ '[' 7925 -lt 10000 ']'
+ FICHIER=velocity_Z_it007925
+ '[' 7925 -lt 1000 ']'
+ '[' 7925 -lt 100 ']'
+ echo velocity_Z_it007925
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007925 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7925+25
+ '[' 7950 -le 8500 ']'
+ '[' 7950 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7950
+ '[' 7950 -lt 100000 ']'
+ FICHIER=velocity_Z_it07950
+ '[' 7950 -lt 10000 ']'
+ FICHIER=velocity_Z_it007950
+ '[' 7950 -lt 1000 ']'
+ '[' 7950 -lt 100 ']'
+ echo velocity_Z_it007950
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007950 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7950+25
+ '[' 7975 -le 8500 ']'
+ '[' 7975 -lt 1000000 ']'
+ FICHIER=velocity_Z_it7975
+ '[' 7975 -lt 100000 ']'
+ FICHIER=velocity_Z_it07975
+ '[' 7975 -lt 10000 ']'
+ FICHIER=velocity_Z_it007975
+ '[' 7975 -lt 1000 ']'
+ '[' 7975 -lt 100 ']'
+ echo velocity_Z_it007975
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it007975 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=7975+25
+ '[' 8000 -le 8500 ']'
+ '[' 8000 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8000
+ '[' 8000 -lt 100000 ']'
+ FICHIER=velocity_Z_it08000
+ '[' 8000 -lt 10000 ']'
+ FICHIER=velocity_Z_it008000
+ '[' 8000 -lt 1000 ']'
+ '[' 8000 -lt 100 ']'
+ echo velocity_Z_it008000
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008000 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8000+25
+ '[' 8025 -le 8500 ']'
+ '[' 8025 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8025
+ '[' 8025 -lt 100000 ']'
+ FICHIER=velocity_Z_it08025
+ '[' 8025 -lt 10000 ']'
+ FICHIER=velocity_Z_it008025
+ '[' 8025 -lt 1000 ']'
+ '[' 8025 -lt 100 ']'
+ echo velocity_Z_it008025
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008025 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8025+25
+ '[' 8050 -le 8500 ']'
+ '[' 8050 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8050
+ '[' 8050 -lt 100000 ']'
+ FICHIER=velocity_Z_it08050
+ '[' 8050 -lt 10000 ']'
+ FICHIER=velocity_Z_it008050
+ '[' 8050 -lt 1000 ']'
+ '[' 8050 -lt 100 ']'
+ echo velocity_Z_it008050
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008050 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8050+25
+ '[' 8075 -le 8500 ']'
+ '[' 8075 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8075
+ '[' 8075 -lt 100000 ']'
+ FICHIER=velocity_Z_it08075
+ '[' 8075 -lt 10000 ']'
+ FICHIER=velocity_Z_it008075
+ '[' 8075 -lt 1000 ']'
+ '[' 8075 -lt 100 ']'
+ echo velocity_Z_it008075
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008075 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8075+25
+ '[' 8100 -le 8500 ']'
+ '[' 8100 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8100
+ '[' 8100 -lt 100000 ']'
+ FICHIER=velocity_Z_it08100
+ '[' 8100 -lt 10000 ']'
+ FICHIER=velocity_Z_it008100
+ '[' 8100 -lt 1000 ']'
+ '[' 8100 -lt 100 ']'
+ echo velocity_Z_it008100
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008100 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8100+25
+ '[' 8125 -le 8500 ']'
+ '[' 8125 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8125
+ '[' 8125 -lt 100000 ']'
+ FICHIER=velocity_Z_it08125
+ '[' 8125 -lt 10000 ']'
+ FICHIER=velocity_Z_it008125
+ '[' 8125 -lt 1000 ']'
+ '[' 8125 -lt 100 ']'
+ echo velocity_Z_it008125
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008125 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8125+25
+ '[' 8150 -le 8500 ']'
+ '[' 8150 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8150
+ '[' 8150 -lt 100000 ']'
+ FICHIER=velocity_Z_it08150
+ '[' 8150 -lt 10000 ']'
+ FICHIER=velocity_Z_it008150
+ '[' 8150 -lt 1000 ']'
+ '[' 8150 -lt 100 ']'
+ echo velocity_Z_it008150
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008150 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8150+25
+ '[' 8175 -le 8500 ']'
+ '[' 8175 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8175
+ '[' 8175 -lt 100000 ']'
+ FICHIER=velocity_Z_it08175
+ '[' 8175 -lt 10000 ']'
+ FICHIER=velocity_Z_it008175
+ '[' 8175 -lt 1000 ']'
+ '[' 8175 -lt 100 ']'
+ echo velocity_Z_it008175
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008175 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8175+25
+ '[' 8200 -le 8500 ']'
+ '[' 8200 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8200
+ '[' 8200 -lt 100000 ']'
+ FICHIER=velocity_Z_it08200
+ '[' 8200 -lt 10000 ']'
+ FICHIER=velocity_Z_it008200
+ '[' 8200 -lt 1000 ']'
+ '[' 8200 -lt 100 ']'
+ echo velocity_Z_it008200
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008200 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8200+25
+ '[' 8225 -le 8500 ']'
+ '[' 8225 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8225
+ '[' 8225 -lt 100000 ']'
+ FICHIER=velocity_Z_it08225
+ '[' 8225 -lt 10000 ']'
+ FICHIER=velocity_Z_it008225
+ '[' 8225 -lt 1000 ']'
+ '[' 8225 -lt 100 ']'
+ echo velocity_Z_it008225
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008225 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8225+25
+ '[' 8250 -le 8500 ']'
+ '[' 8250 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8250
+ '[' 8250 -lt 100000 ']'
+ FICHIER=velocity_Z_it08250
+ '[' 8250 -lt 10000 ']'
+ FICHIER=velocity_Z_it008250
+ '[' 8250 -lt 1000 ']'
+ '[' 8250 -lt 100 ']'
+ echo velocity_Z_it008250
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008250 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8250+25
+ '[' 8275 -le 8500 ']'
+ '[' 8275 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8275
+ '[' 8275 -lt 100000 ']'
+ FICHIER=velocity_Z_it08275
+ '[' 8275 -lt 10000 ']'
+ FICHIER=velocity_Z_it008275
+ '[' 8275 -lt 1000 ']'
+ '[' 8275 -lt 100 ']'
+ echo velocity_Z_it008275
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008275 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8275+25
+ '[' 8300 -le 8500 ']'
+ '[' 8300 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8300
+ '[' 8300 -lt 100000 ']'
+ FICHIER=velocity_Z_it08300
+ '[' 8300 -lt 10000 ']'
+ FICHIER=velocity_Z_it008300
+ '[' 8300 -lt 1000 ']'
+ '[' 8300 -lt 100 ']'
+ echo velocity_Z_it008300
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008300 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8300+25
+ '[' 8325 -le 8500 ']'
+ '[' 8325 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8325
+ '[' 8325 -lt 100000 ']'
+ FICHIER=velocity_Z_it08325
+ '[' 8325 -lt 10000 ']'
+ FICHIER=velocity_Z_it008325
+ '[' 8325 -lt 1000 ']'
+ '[' 8325 -lt 100 ']'
+ echo velocity_Z_it008325
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008325 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8325+25
+ '[' 8350 -le 8500 ']'
+ '[' 8350 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8350
+ '[' 8350 -lt 100000 ']'
+ FICHIER=velocity_Z_it08350
+ '[' 8350 -lt 10000 ']'
+ FICHIER=velocity_Z_it008350
+ '[' 8350 -lt 1000 ']'
+ '[' 8350 -lt 100 ']'
+ echo velocity_Z_it008350
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008350 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8350+25
+ '[' 8375 -le 8500 ']'
+ '[' 8375 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8375
+ '[' 8375 -lt 100000 ']'
+ FICHIER=velocity_Z_it08375
+ '[' 8375 -lt 10000 ']'
+ FICHIER=velocity_Z_it008375
+ '[' 8375 -lt 1000 ']'
+ '[' 8375 -lt 100 ']'
+ echo velocity_Z_it008375
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008375 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8375+25
+ '[' 8400 -le 8500 ']'
+ '[' 8400 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8400
+ '[' 8400 -lt 100000 ']'
+ FICHIER=velocity_Z_it08400
+ '[' 8400 -lt 10000 ']'
+ FICHIER=velocity_Z_it008400
+ '[' 8400 -lt 1000 ']'
+ '[' 8400 -lt 100 ']'
+ echo velocity_Z_it008400
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008400 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8400+25
+ '[' 8425 -le 8500 ']'
+ '[' 8425 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8425
+ '[' 8425 -lt 100000 ']'
+ FICHIER=velocity_Z_it08425
+ '[' 8425 -lt 10000 ']'
+ FICHIER=velocity_Z_it008425
+ '[' 8425 -lt 1000 ']'
+ '[' 8425 -lt 100 ']'
+ echo velocity_Z_it008425
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008425 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8425+25
+ '[' 8450 -le 8500 ']'
+ '[' 8450 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8450
+ '[' 8450 -lt 100000 ']'
+ FICHIER=velocity_Z_it08450
+ '[' 8450 -lt 10000 ']'
+ FICHIER=velocity_Z_it008450
+ '[' 8450 -lt 1000 ']'
+ '[' 8450 -lt 100 ']'
+ echo velocity_Z_it008450
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008450 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8450+25
+ '[' 8475 -le 8500 ']'
+ '[' 8475 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8475
+ '[' 8475 -lt 100000 ']'
+ FICHIER=velocity_Z_it08475
+ '[' 8475 -lt 10000 ']'
+ FICHIER=velocity_Z_it008475
+ '[' 8475 -lt 1000 ']'
+ '[' 8475 -lt 100 ']'
+ echo velocity_Z_it008475
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008475 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8475+25
+ '[' 8500 -le 8500 ']'
+ '[' 8500 -lt 1000000 ']'
+ FICHIER=velocity_Z_it8500
+ '[' 8500 -lt 100000 ']'
+ FICHIER=velocity_Z_it08500
+ '[' 8500 -lt 10000 ']'
+ FICHIER=velocity_Z_it008500
+ '[' 8500 -lt 1000 ']'
+ '[' 8500 -lt 100 ']'
+ echo velocity_Z_it008500
+ /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/bin/xcombine_vol_data 0 31 velocity_Z_it008500 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/OUTPUT_FILES/DATABASES_MPI /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie 0
error opening Par_file
+ it=8500+25
+ '[' 8525 -le 8500 ']'
+ tar -jcvf /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie.tar.bz2 /ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/DSM_FOR_SPECFEM3D/EXAMPLES/example_simple_small/movie
tar: Removing leading `/' from member names
+ cd ..
++ date
+ echo Wed Jul 30 23:03:01 CEST 2014
+ exit 0
