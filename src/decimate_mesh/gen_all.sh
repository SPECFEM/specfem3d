#!/bin/sh
# Generates the mesh, nodes_coords and mat files from a given Abaqus one. 
# The "filename" variable holds the name of the Abaqus file.

filename="HOMO_3D_lisse_300"
nbElems=$(tail -1 $filename | awk -F\, '{print $1}')  
echo $nbElems > mesh
nb=$(cat mesh)
tail -$nb $filename | awk -F\, '{print $2 $3 $4 $5 $6 $7 $8 $9}' >> mesh
nbLines=$(wc -l $filename | awk '{print $1}') 
nbNodes=$(($nbLines-$nbElems-4))
echo $nbNodes > nodes_coords
head -$(($nbNodes+3)) $filename | tail -$nbNodes |  awk -F\, '{print $2 $3 $4}' >> nodes_coords
echo 1 > mat
for i in `seq 2 $nbNodes`;
do
    echo 1 >> mat
done    
