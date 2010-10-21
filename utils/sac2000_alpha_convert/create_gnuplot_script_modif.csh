#!/bin/csh

echo "set term x11"
echo "#set xrange [0:1400]"

foreach file ($*)

set newfile = `basename $file .DW.BHE.dataa_ascii.sac.veloc.ascii`
echo "plot '$file' w l 1, '$newfile.DW.BHN.dataa_ascii.sac.veloc.ascii' w l 3, '$newfile.DW.BHZ.dataa_ascii.sac.veloc.ascii' w l 4"
echo 'pause -1 "hit key"'

end

