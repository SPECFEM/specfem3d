#!/bin/csh

# compile and run the mesher and then the solver, then collect the seismograms

set OUTDIR="SEM"

sleep 1
make clean
sleep 1
make meshfem3D
sleep 10
go_mesher

sleep 5
make clean
sleep 1
make specfem3D
sleep 10
go_solver

sleep 5
if (! -d $OUTDIR) then
   mkdir $OUTDIR
endif 
cd $OUTDIR
cp ../DATA/CMTSOLUTION .
cp ../DATA/STATIONS .
cp ../DATA/Par_file .
collect_seismos < ../collect_seismos.in

