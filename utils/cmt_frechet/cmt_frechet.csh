#!/bin/csh
#
# NOTE: this script uses a local script "collect_seismos" used at Caltech.
#             Please modify for your own purposes.
#
sleep 1
make clean
sleep 1
make generate_databases
sleep 1
make specfem3D

sleep 10
go_generate_databases

# save the CMTSOLUTION and the STATIONS file
mkdir -p CMT_FRECHET
cp in_data_files/CMTSOLUTION CMT_FRECHET
cp in_data_files/STATIONS CMT_FRECHET

# calculate synthetics
sleep 180
go_solver_cmt

# collect the synthetics
cd CMT_FRECHET
collect_seismos < ../collect_seismos.in
foreach file (*.semd)
  set shortfile = `basename $file .semd `
  mv $file ${shortfile}
end
cd ..

# make the CMTSOLUTION files needed for the calculation of Frechet derivatives
cd in_data_files
# backup the CMTSOLUTION file
cp CMTSOLUTION CMTSOLUTION_BACKUP
./xmake_cmtsolution_files
cd ..

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

