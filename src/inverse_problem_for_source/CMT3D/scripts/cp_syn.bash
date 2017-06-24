#!/bin/bash
# this is an old script used to copy and untar derivative syns,
# copy measurement files, add up syns, change file extension,
# copy inversion pars, figure out the dlat, dlon from CMTs.
if [ $# != 1 ]; then
  echo "Usage: cp_syn.bash evid-file"; exit
fi

extra_dir=""

for evid in `cat $1 `; do
  echo "==== $evid ===="
  mkdir -p ${evid}${extra_dir};
  cd ${evid}${extra_dir}; mkdir -p syn; mkdir -p backup
## copy sem syn from cluster
#  scp $cluster/${evid}${extra_dir}/* syn/
  mv ../syns/*${evid}*.bz syn/
  cp ../event_info/CMTSOLUTION_$evid CMTSOLUTION; cp CMTSOLUTION syn/
  cp ../event_info/STATIONS .
## Unzip synthetics
  cd syn
  tar xjvf CMTSOLUTION_${evid}.tar.bz
  for ext in Mrr Mtt Mpp Mrt Mrp Mtp dep lat lon; do
    echo " Untar $ext ..."; tar xjvf ${evid}_${ext}.tar.bz >/dev/null
    if [ $? != 0 ]; then
      echo "Error untaring ... $ext ..."; exit
    fi
    mv -f ${evid}_${ext}.tar.bz ../backup
    rename "s/semd.sac.$ext/$ext/" *.semd.sac.$ext
  done
  ls -1 *.Mrr | perl -pi -e 's/\.Mrr//g' -- > file_list
  xadd_frechet_derivatives s file_list ../CMTSOLUTION 1.0e22
  cd ..
## Copy measurement file
  cp ../event_info/MEASUREMENT_WINDOWS_${evid}_T00?_T030_m12 .
## here check if MEASURE file has the right number of data/syn pairsa
  mafile="MEASUREMENT_WINDOWS_${evid}_ALL"
  nma=0
  rm -f $mafile
  for freq in 2 3 6; do
    mfile="MEASUREMENT_WINDOWS_${evid}_T00${freq}_T030_m12"
    n0=`grep "^[DS]" $mfile | wc | awk '{print $1}'`
    n1=`echo "$n0 / 2" | bc`
    n2=`head -1 $mfile`
    if [ $n1 != $n2 ]; then
      echo "Number of files inconsistent in $mfile: $n1 and $n2"; exit
    fi
## here modify dir names: DATA -> data_T006_T030, SYN -> syn_T006_T030
    perl -pi -e "s/DATA/data_T00${freq}_T030/g" $mfile
    perl -pi -e "s/SYN/syn_T00${freq}_T030/g" $mfile
    perl -pi -e "s/\.T00${freq}_T030//g" $mfile
    perl -pi -e "s/\.semd\.sac\.m12//g" $mfile
    awk 'NR > 1 {print $0}' $mfile >> $mafile
    nm0=`head -n 1 $mfile`
    nma=`echo $nm0 + $nma | bc `
  done
## right number of measurements
  echo $nma | cat - $mafile > out.tmp; mv out.tmp $mafile
#### Link cmt3d and grid3d_flexwin
  ln -s ../../grd_cmt3d/cmt3d/cmt3d_flexwin
  ln -s ../../grd_cmt3d/grid3d/grid3d_flexwin
## here modify INVERSION.PAR GRID3D.PAR for MEASURE file
  cp ../../grd_cmt3d/cmt3d/INVERSION.PAR .
  cp ../../grd_cmt3d/grid3d/GRID3D.PAR .
  perl -pi -e "s/flexwin\.out/$mafile/g" INVERSION.PAR
  mv INVERSION.PAR INVERSION.PAR.SAVE
  perl -pi -e "s/flexwin\.out/$mafile/g" GRID3D.PAR
##  make sure you have the right dmoment, ddepth, and dlocation
  dlat=`diff -b syn/CMTSOLUTION  syn/CMTSOLUTION_lat  | grep lat | awk '{print $3}' | perl -pi -e 's/\n/ /g' | awk '{print $2-$1}'`
  ddep=`diff -b syn/CMTSOLUTION  syn/CMTSOLUTION_dep  | grep dep | awk '{print $3}' | perl -pi -e 's/\n/ /g' | awk '{print $2-$1}'`
  perl -pi -e "s/^.*1.0e22*/$dlat $ddep 1.0e22/" INVERSION.PAR.SAVE
####
  cd ..;
done
