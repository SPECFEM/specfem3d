#?/bin/bash

numsrc=$1

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

cd Source_${numsrc}
cp ../*.x .
cp ../submit_reconst2.sh .
cp ../reconstruction2.par reconstruction.par
cp ../input_box.txt .

cp input_box.txt MXX_P_MYY/.
cp input_box.txt MXY_MXX_M_MYY/.
cp input_box.txt MXZ_MYZ/.
cp input_box.txt MZZ/.

sbatch submit_reconst2.sh

cd ..


