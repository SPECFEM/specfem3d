#!/bin/bash -f

numrec=$1
colat=0.
incr=`echo $numrec | awk '{print (180./( $1 - 1 ))}'`

echo $numrec > receivers_$1.dat
ind=0

echo $numrec $colat $incr

while [  $ind -lt $numrec ]; do

let ind=$ind+1
echo $colat $ind
echo $colat "0.0" >> receivers_$1.dat
colat=`echo $colat $incr | awk '{print ($1 + $2) }'`

done