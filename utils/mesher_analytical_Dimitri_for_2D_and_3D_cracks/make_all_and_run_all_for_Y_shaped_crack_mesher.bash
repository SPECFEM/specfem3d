#!/bin/bash

rm -f xtutu1 xtutu2 xtutu3 *__genmod.*

ifort -std08 -implicitnone -warn all -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv -o xtutu1 mesh_a_Y_shaped_crack_oriented_upwards.f90

ifort -std08 -implicitnone -warn all -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv -o xtutu2 create_database_files_for_external_faces_of_the_3D_mesh.f90

ifort -std08 -implicitnone -warn all -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv -o xtutu3 create_OpenDX_files_to_view_the_3D_mesh.f90

rm -f *__genmod.*

mkdir -p DATA
mkdir -p DATA_3D

# run the two codes to create the 2D mesh, the 3D mesh, and the database files for SPECFEM3D_Cartesian
./xtutu1
echo "creating the database files for the outer faces of the 3D..."
./xtutu2
echo "Done."

echo " "
echo "Now creating an OpenDX file to check the 3D mesh created..."
./xtutu3
echo "Done."

