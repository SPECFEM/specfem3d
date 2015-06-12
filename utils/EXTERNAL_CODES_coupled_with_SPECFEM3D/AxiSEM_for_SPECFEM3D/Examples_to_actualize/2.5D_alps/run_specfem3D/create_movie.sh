#!/bin/bash
source ./set_path_and_params.sh 
IN=./OUTPUT_FILES/DATABASES_MPI/
OUT=./OUTPUT_FILES/VOLUME_VIDEO

mkdir ${OUT}
#cd bin

NAME_FIELD="velocity_Z_it"

declare -i it itmax istep nproc                                                                                                                              
ires=1 # 1 high resolution moive, 0 low relsolution

itmax=`grep NSTEP DATA/Par_file | cut -c 35-` # NSTEP
istep=`grep NTSTEP_BETWEEN_FRAMES DATA/Par_file | cut -c 35-` # NTSTEP_BETWEEN_FRAMES 

nproc=$NPROC
nproc="$nproc-1"

#itmax=6000
#istep=50
#nproc=15
it=${istep};
 it=800                                                
while [ "$it" -le "$itmax" ] ; do                                                                                                   
echo "it=" ${it}
  if [[ "$it" -lt 10 ]]; then                                                                                                       
    TAG="00000"${it}       
  fi;                                                                                                                               

  if [[ "$it" -ge 10 && "$it" -lt 100 ]]; then                                                                                      
    TAG="0000"${it}  
  fi;                                                                                                                               

  if [[ "$it" -ge 100 && "$it" -lt 1000 ]]; then                                                                                    
    TAG="000"${it}     
  fi;                                                                                                                               

  if [[ "$it" -ge 1000 && "$it" -lt 10000 ]]; then                                                                                  
    TAG="00"${it}
  fi;                                                                                                                               

  if [[ "$it" -ge 10000 && "$it" -lt 100000 ]]; then
    TAG="0"${it}
  fi;

  fichier=$NAME_FIELD$TAG
  echo $fichier
  $SEMBIN/xcombine_vol_data 0 ${nproc} $fichier $IN $OUT $ires
  #mesh2vtu.pl -i $OUT/$fichier".mesh" -o $OUT/$fichier".vtu"
  it="$it+$istep"
 
done;



