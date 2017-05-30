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
cp DATA/CMTSOLUTION CMT_FRECHET
cp DATA/STATIONS CMT_FRECHET

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
cd DATA
# backup the CMTSOLUTION file
cp CMTSOLUTION CMTSOLUTION_BACKUP
./xmake_cmtsolution_files
cd ..

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

