
#set term postscript color solid
#set output "results_of_seismogram_comparison.ps"

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'REF_SEIS/DB.X20.BXZ.semd' w l lc 1, 'OUTPUT_FILES/DB.X20.BXZ.semd' w l lc 3
pause -1 "hit key..."

plot 'REF_SEIS/DB.X30.BXZ.semd' w l lc 1, 'OUTPUT_FILES/DB.X30.BXZ.semd' w l lc 3
pause -1 "hit key..."

plot 'REF_SEIS/DB.X40.BXZ.semd' w l lc 1, 'OUTPUT_FILES/DB.X40.BXZ.semd' w l lc 3
pause -1 "hit key..."

plot 'REF_SEIS/DB.X50.BXZ.semd' w l lc 1, 'OUTPUT_FILES/DB.X50.BXZ.semd' w l lc 3
pause -1 "hit key..."

