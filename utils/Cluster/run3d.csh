#!/bin/csh

# compile and run the mesher and then the solver, then collect the seismograms

set OUTDIR="in_out_files/SEM"

sleep 1
make clean
sleep 1
make xgenerate_databases
sleep 10
go_generate_databases

sleep 5
make clean
sleep 1
make xspecfem3D
sleep 10
go_solver

sleep 5
if (! -d $OUTDIR) then
   mkdir $OUTDIR
endif
cd $OUTDIR
cp ../../in_data_files/CMTSOLUTION .
cp ../../in_data_files/STATIONS .
cp ../../in_data_files/Par_file .
collect_seismos < ../collect_seismos.in

