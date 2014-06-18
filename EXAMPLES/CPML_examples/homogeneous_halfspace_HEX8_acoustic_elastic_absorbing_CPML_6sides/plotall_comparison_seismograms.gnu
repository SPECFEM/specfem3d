
#set term postscript color solid
#set output "results_of_seismogram_comparison.ps"

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'REF_SEIS/X20.DB.BXZ.semd' w l lc 1, 'OUTPUT_FILES/X20.DB.BXZ.semd' w l lc 3
pause -1 "hit key..."

plot 'REF_SEIS/X30.DB.BXZ.semd' w l lc 1, 'OUTPUT_FILES/X30.DB.BXZ.semd' w l lc 3
pause -1 "hit key..."

plot 'REF_SEIS/X40.DB.BXZ.semd' w l lc 1, 'OUTPUT_FILES/X40.DB.BXZ.semd' w l lc 3
pause -1 "hit key..."

plot 'REF_SEIS/X50.DB.BXZ.semd' w l lc 1, 'OUTPUT_FILES/X50.DB.BXZ.semd' w l lc 3
pause -1 "hit key..."

