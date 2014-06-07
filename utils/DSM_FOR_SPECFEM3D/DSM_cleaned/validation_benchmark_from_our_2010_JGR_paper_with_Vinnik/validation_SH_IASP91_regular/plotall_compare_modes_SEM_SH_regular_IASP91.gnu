set term wxt
set term x11

#set term postscript eps monochrome dashed "Helvetica" 22
set term postscript eps color solid "Helvetica" 22
set output "compare_modes_SEM_SH_regular_IASP91.eps"

set xlabel "Time (s)"
set ylabel "Displacement (m)"

set xrange [1350:1750]

set label 1 at 1530,+0.0012
set label 1 "SHdiff"

plot 'S100_00.DK.LHN.semd.convolved_regular_IASP91_filtered.dat' t 'SH SEM IASP91' w l 1, 'DSM_SH_regular_filtered.dat.convolved' t 'SH DSM IASP91' w l 2, 'seismogram_Gemini_LHE.txt_iasp91_regular_nocrust_noQ_smoothed_filtered.convolved' t 'SH Gemini IASP91' w l 3, 'Vinnik.S100.ASC_comp3.convolved_rescaled_regular_IASP91_filtered.dat' t 'SH Modes IASP91' w l 4
pause -1 "hit key"

