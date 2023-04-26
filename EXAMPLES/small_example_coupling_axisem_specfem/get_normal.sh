#!/bin/bash

# number of processors
declare -i NPROC
NPROC=$1
NPROC="$NPROC-1"

# database directory
DATABASES_MPI=DATABASES_MPI

# create output file
outfile="normals.txt"
rm $outfile
touch $outfile

for I in `seq 0 1 $NPROC`;
do
 i=$(printf "%06d" $I)
 echo "READING :" $DATABASES_MPI"/proc"$i"_normal.txt"
 infile=$DATABASES_MPI"/proc"$i"_normal.txt"
 awk  -v II=$I ' NF>5 {print $0  "   " II}' < $infile  >> $outfile
done

echo
echo " Wrote  " $outfile " file "
echo
echo " X, Y, Z, NX, NY, NZ, iproc"
echo
echo
