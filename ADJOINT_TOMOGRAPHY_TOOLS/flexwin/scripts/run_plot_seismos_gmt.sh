#!/bin/sh

dir=$1

for win in ${dir}/*.win ; do
  echo $win
  base=`echo $win | sed s/.win//`
  ./plot_seismos_stalta.gmt $base
done

./make_single_pdf.sh $dir ${dir}.pdf
