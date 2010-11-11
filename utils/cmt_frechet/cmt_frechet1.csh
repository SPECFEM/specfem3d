#!/bin/csh

# set up frechet derivative calculations for the CMTSOLUTION
# Qinya Liu, May 2007, Caltech
#
# NOTE: this script uses a local script "collect_seismos" used at Caltech.
#             Please modify for your own purposes.
#

foreach extension (Mrr Mtt Mpp Mrt Mrp Mtp depth latitude longitude)

# setup the calculation for the Frechet derivative with respect to ${extension}
rm in_data_files/CMTSOLUTION
mv in_data_files/CMTSOLUTION_${extension} in_data_files/CMTSOLUTION

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
rm in_data_files/CMTSOLUTION
mv in_data_files/CMTSOLUTION_BACKUP in_data_files/CMTSOLUTION

