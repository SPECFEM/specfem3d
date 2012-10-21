
set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'REF_SEIS/X1.DB.BXZ.semd' w l lc 1, 'in_out_files/OUTPUT_FILES/X1.DB.BXZ.semd' w l lc 3
pause -1 "hit key..."

plot 'REF_SEIS/X10.DB.BXZ.semd' w l lc 1, 'in_out_files/OUTPUT_FILES/X10.DB.BXZ.semd' w l lc 3
pause -1 "hit key..."

plot 'REF_SEIS/X20.DB.BXZ.semd' w l lc 1, 'in_out_files/OUTPUT_FILES/X20.DB.BXZ.semd' w l lc 3
pause -1 "hit key..."

plot 'REF_SEIS/X30.DB.BXZ.semd' w l lc 1, 'in_out_files/OUTPUT_FILES/X30.DB.BXZ.semd' w l lc 3
pause -1 "hit key..."

plot 'REF_SEIS/X40.DB.BXZ.semd' w l lc 1, 'in_out_files/OUTPUT_FILES/X40.DB.BXZ.semd' w l lc 3
pause -1 "hit key..."

plot 'REF_SEIS/X50.DB.BXZ.semd' w l lc 1, 'in_out_files/OUTPUT_FILES/X50.DB.BXZ.semd' w l lc 3
pause -1 "hit key..."

