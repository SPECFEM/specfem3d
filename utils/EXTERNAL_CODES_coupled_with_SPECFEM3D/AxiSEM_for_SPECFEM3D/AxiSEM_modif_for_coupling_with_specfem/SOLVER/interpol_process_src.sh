#?/bin/bash

module load intel/15.0.0.090
module load bullxmpi/1.2.8.3


mode=$2
numsrc=$1
it=$3
dt=0.02
mytype=0

tshift=`grep "time shift:"  Source_${numsrc}/CMTSOLUTION  |awk '{print $3}'`
hdur=`grep "half duration:" Source_${numsrc}/CMTSOLUTION  |awk '{print $3}'`

rm -f output input fsource.bin
rm -f Source_${numsrc}/fsource.bin
rm -f Source_${numsrc}/verif.bin
rm -f Source_${numsrc}/verifconv.bin
rm -f Source_${numsrc}/stfint.bin
rm -f Source_${numsrc}/stalta.bin

if [ $mode == '1' ]; then

  echo "$mytype $hdur $tshift $dt" > input
  ./source_lithos.x < input > output

  nt=`cat output | grep -v "type,hdur,shift,dt" | awk '{print $1}'`

        lat=`grep "latitude:"  Source_${numsrc}/CMTSOLUTION  |awk '{print $2}'`
        lon=`grep "longitude:" Source_${numsrc}/CMTSOLUTION  |awk '{print $2}'`

        echo "SEM"            > reconstruction2.par
        echo "input_box.txt" >> reconstruction2.par
        echo "480"           >> reconstruction2.par
        echo "$lat $lon"     >> reconstruction2.par
        echo "44.6 7 30"     >> reconstruction2.par
        echo "4"             >> reconstruction2.par

        sed "s|mondossier|Source_$numsrc|"  submit_reconst.sh >  submit_reconst2.sh
        sed -i "s|monnom|$numsrc/RECON|"    submit_reconst2.sh

  echo "Estimate starting time and compute stf"
  echo "SEM"                            > interpolation2.par
  echo "../../../../Goodone" >> interpolation2.par
  echo "20 400 1.5"                      >> interpolation2.par
  echo "0 400"                         >> interpolation2.par
  echo "1 $nt $dt"                   >> interpolation2.par
  echo "fsource.bin"                   >> interpolation2.par
  echo "480"                           >> interpolation2.par
  echo "0.2"                           >> interpolation2.par


  cp fsource.bin Source_${numsrc}/.

        cd Source_${numsrc}

        cp ../*.x .
        cp ../interpolation2.par interpolation.par
        cp ../reconstruction2.par reconstruction.par
        cp ../submit_eval.sh .
  cp ../input .
  cp ../output .

        sed "s|mondossier|Source_$numsrc|"  submit_eval.sh >  submit_eval2.sh
        sed -i "s|monnom|$numsrc/INTERP|"    submit_eval2.sh

        sbatch submit_eval2.sh

        cd ..


elif [ $mode == '2' ]; then

        echo "Initerpolate source"
    echo "$mytype $hdur $tshift $dt" > input
        ./source_lithos.x < input > output
  nt=`cat output | grep -v "type,hdur,shift,dt" | awk '{print $1}'`

  cp fsource.bin Source_${numsrc}/.

  lat=`grep "latitude:"  Source_${numsrc}/CMTSOLUTION  |awk '{print $2}'`
  lon=`grep "longitude:" Source_${numsrc}/CMTSOLUTION  |awk '{print $2}'`

  echo "SEM"            > reconstruction2.par
  echo "input_box.txt" >> reconstruction2.par
  echo "480"           >> reconstruction2.par
  echo "$lat $lon"     >> reconstruction2.par
  echo "44.6 7 30"     >> reconstruction2.par
  echo "4"             >> reconstruction2.par

  sed "s|mondossier|Source_$numsrc|"  submit_reconst.sh >  submit_reconst2.sh
  sed -i "s|monnom|$numsrc/RECON|"    submit_reconst2.sh

        echo "SEM"                            > interpolation2.par
        echo "../../../../Goodone" >> interpolation2.par
        echo "20 400 1.5"                      >> interpolation2.par
        echo "1 $3"                          >> interpolation2.par
        echo "1 $nt $dt"                   >> interpolation2.par
        echo "fsource.bin"                   >> interpolation2.par
        echo "480"                           >> interpolation2.par
        echo "0.2"                           >> interpolation2.par

  sed "s|mondossier|Source_$numsrc|"  submit_interp.sh > submit_interp2.sh
  sed -i "s|monnom|$numsrc/INTERP|"    submit_interp2.sh

  cd Source_${numsrc}
  cp ../*.x .
  cp ../submit_interp2.sh .
  cp ../interpolation2.par interpolation.par
  cp ../reconstruction2.par reconstruction.par
  cp ../input .
  cp ../output .
  sbatch submit_interp2.sh

  cd ..
fi


