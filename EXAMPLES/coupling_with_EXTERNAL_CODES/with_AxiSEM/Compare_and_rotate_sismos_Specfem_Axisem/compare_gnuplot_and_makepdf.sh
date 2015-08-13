#!/bin/bash

### Launch the script from the directory containing subdirs SISMOS_RUN_1 and SISMOS_RUN_2

mkdir SISMOS_PLOT_COMPARED_Axisem_and_Specfem_converted

gnuplot allgn.plt

mv *.eps SISMOS_PLOT_COMPARED_Axisem_and_Specfem_converted
cd SISMOS_PLOT_COMPARED_Axisem_and_Specfem_converted

epstopdf SY.S0001.E.semd_COMPARED.eps
epstopdf SY.S0001.N.semd_COMPARED.eps
epstopdf SY.S0001.Z.semd_COMPARED.eps
epstopdf SY.S0127.E.semd_COMPARED.eps
epstopdf SY.S0127.N.semd_COMPARED.eps
epstopdf SY.S0127.Z.semd_COMPARED.eps
epstopdf SY.S0230.E.semd_COMPARED.eps
epstopdf SY.S0230.N.semd_COMPARED.eps
epstopdf SY.S0230.Z.semd_COMPARED.eps
epstopdf SY.S0253.E.semd_COMPARED.eps
epstopdf SY.S0253.N.semd_COMPARED.eps
epstopdf SY.S0253.Z.semd_COMPARED.eps
epstopdf SY.S0076.E.semd_COMPARED.eps
epstopdf SY.S0076.N.semd_COMPARED.eps
epstopdf SY.S0076.Z.semd_COMPARED.eps

rm *.eps

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=concatall.pdf *.pdf

