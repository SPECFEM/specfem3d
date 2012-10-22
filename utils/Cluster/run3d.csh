#!/bin/csh

# compile and run the database generator and then the solver

set OUTDIR="OUTPUT_FILES"

sleep 1
make clean
sleep 1
make xgenerate_databases
sleep 1
make xspecfem3D

sleep 10
go_generate_databases
sleep 1
go_solver

sleep 5
if (! -d $OUTDIR) then
   mkdir $OUTDIR
endif
cd $OUTDIR
cp ../../DATA/CMTSOLUTION .
cp ../../DATA/STATIONS .
cp ../../DATA/Par_file .

