#!/bin/csh 

   sed -e '1,$s/D-00/e-/g' < gmt_shaking_hollywood_veloc.irreg > tutu
   mv tutu gmt_shaking_hollywood_veloc.irreg
   sed -e '1,$s/  / /g' < gmt_shaking_hollywood_veloc.irreg > tutu
   mv tutu gmt_shaking_hollywood_veloc.irreg
   sed -e '1,$s/  / /g' < gmt_shaking_hollywood_veloc.irreg > tutu
   mv tutu gmt_shaking_hollywood_veloc.irreg
   sed -e '1,$s/  / /g' < gmt_shaking_hollywood_veloc.irreg > tutu
   mv tutu gmt_shaking_hollywood_veloc.irreg

   sed -e '1,$s/D-00/e-/g' < gmt_shaking_hollywood_displ.irreg > tutu
   mv tutu gmt_shaking_hollywood_displ.irreg
   sed -e '1,$s/  / /g' < gmt_shaking_hollywood_displ.irreg > tutu
   mv tutu gmt_shaking_hollywood_displ.irreg
   sed -e '1,$s/  / /g' < gmt_shaking_hollywood_displ.irreg > tutu
   mv tutu gmt_shaking_hollywood_displ.irreg
   sed -e '1,$s/  / /g' < gmt_shaking_hollywood_displ.irreg > tutu
   mv tutu gmt_shaking_hollywood_displ.irreg

