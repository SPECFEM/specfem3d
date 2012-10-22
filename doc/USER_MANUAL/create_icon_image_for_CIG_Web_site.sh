#!/bin/sh

# Dimitri Komatitsch, November 2010: use "convert" from the Imagemagick package

convert -quality 100 -resize 38% -colorspace RGB figures/specfem_3d_Cartesian-cover.pdf jpg:specfem_3d_Cartesian-cover_small.jpg

