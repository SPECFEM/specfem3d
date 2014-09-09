#!/bin/sh

dir=$1
out_pdf=$2

pdftk ${dir}/*stalta.pdf cat output $out_pdf
pdfnup --nup 3x3 $out_pdf

