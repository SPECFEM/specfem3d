#!/bin/csh

# code to compute attenuation constants for SPECFEM3D_SESAME code

rm -f xcompute
gcc -o xcompute attenuation_sesame.c -lm

echo "code has been compiled"
echo "run it with:  xcompute "

