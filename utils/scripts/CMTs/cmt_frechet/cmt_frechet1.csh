#!/bin/csh

# set up frechet derivative calculations for the CMTSOLUTION
# Qinya Liu, May 2007, Caltech
#
# NOTE: this script uses a local script "collect_seismos" used at Caltech.
#             Please modify for your own purposes.
#

foreach extension (Mrr Mtt Mpp Mrt Mrp Mtp depth latitude longitude)

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

