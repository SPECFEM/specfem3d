set term wxt
set term x11 

#set term postscript eps monochrome dashed "Helvetica" 22
set term postscript eps color solid "Helvetica" 22

set xlabel "Time (s)"
set ylabel "Displacement (m)"

###################################################

set output "compare_modes_SEM_SV_regular_IASP91_zoom.eps"

set xrange [1420:1590]

set yrange [-0.00335:+0.0024]

set label 1 at 1523.3,+0.0007
set label 1 "SVdiff"

set label 2 at 1467.2,-0.0023
set label 2 "SKS"

set label 3 at 1487,+0.0011
set label 3 "SKKS"

plot 'S100_00.DK.LHE.semd.convolved_resent_by_Vinnik_filtered' t 'SV SEM IASP91' w l 1, 'DSM_SV_regular_filtered.dat.convolved' t 'SV DSM IASP91' w l 2, 'seismogram_Gemini_LHN.txt_iasp91_regular_nocrust_noQ_smoothed_filtered.convolved' t 'SV Gemini IASP91' w l 3
pause -1 "hit key"

############################################################

set output "compare_modes_SEM_SV_regular_IASP91.eps"

set xrange [1350:1750]

set yrange [-0.012:+0.017]

set label 1 at 1500,+0.0025
set label 1 "SVdiff"

set label 2 at 1610,+0.007
set label 2 "PS"

set label 3 at 1650,-0.008
set label 3 "PPS"

set label 4 at 1448,-0.0045
set label 4 "SKS"

plot 'S100_00.DK.LHE.semd.convolved_resent_by_Vinnik_filtered' t 'SV SEM IASP91' w l 1, 'DSM_SV_regular_filtered.dat.convolved' t 'SV DSM IASP91' w l 2, 'seismogram_Gemini_LHN.txt_iasp91_regular_nocrust_noQ_smoothed_filtered.convolved' t 'SV Gemini IASP91' w l 3
pause -1 "hit key"

############################################################


