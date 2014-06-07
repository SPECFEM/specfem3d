PREFIX=velocity_Z_it
OUT=./movie


declare -i it itmax istep


istep=25
itmax=3500 #30000


it=0

while [ "$it" -le "$itmax" ] ; do

   if [ "$it" -lt 1000000 ]; then 
     FICHIER=$PREFIX${it}
   fi;

   if [ "$it" -lt 100000 ]; then
    FICHIER=$PREFIX"0"${it}
   fi;

   if [ "$it" -lt 10000 ]; then
      FICHIER=$PREFIX"00"${it}
   fi;

   if [ "$it" -lt 1000 ]; then
     FICHIER=$PREFIX"000"${it}
   fi;

   if [ "$it" -lt 100 ]; then
      FICHIER=$PREFIX"0000"${it}
   fi;


   echo $FICHIER.mesh
   ./mesh2vtu.pl -i $OUT/$FICHIER.mesh -o $FICHIER.vtu
   it="$it+$istep"

done;

echo $PREFIX

mv $PREFIX* ./movie_vtu/


