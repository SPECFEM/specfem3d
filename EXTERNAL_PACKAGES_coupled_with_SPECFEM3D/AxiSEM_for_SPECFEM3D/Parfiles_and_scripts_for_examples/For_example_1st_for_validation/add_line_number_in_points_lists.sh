#!/bin/bash

function add_numb_pts_and_delete_if_exist ()
{
if [ -e $1 ] ; then
  NUMBOFPOINTS=$(wc -l < $1)
  echo ${NUMBOFPOINTS} | cat - $1 > $2
  rm -f $1
else
  echo "Files ""list_ggl_boundary..."" have not been copied or are already deleted ==> Check files ""input_box..."", they may be already OK !!"
fi
}

add_numb_pts_and_delete_if_exist list_ggl_boundary_spherical.txt input_box.txt

add_numb_pts_and_delete_if_exist list_ggl_boundary_cartesian.txt input_box_sem_cart.txt

