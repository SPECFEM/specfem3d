#!/bin/csh

foreach extension (Mrt Mrp Mtp)

# setup the calculation for the Frechet derivative with respect to ${extension}
rm DATA/CMTSOLUTION
mv DATA/CMTSOLUTION_${extension} DATA/CMTSOLUTION

sleep 180
go_solver_cmt

cd CMT_FRECHET
collect_seismos < ../collect_seismos.in
foreach file (*.semd)
  set shortfile = `basename $file .semd `
  mv $file ${shortfile}.${extension}
end
cd ..
 
end

# put the backup CMTSOLUTION file back
rm DATA/CMTSOLUTION
mv DATA/CMTSOLUTION_BACKUP DATA/CMTSOLUTION

