#!/usr/local/bin/gnuplot -persist
#
#    

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'X1.DB.BXZ.semd' w l,'X10.DB.BXZ.semd' w l,'X20.DB.BXZ.semd' w l,'X30.DB.BXZ.semd' w l,'X40.DB.BXZ.semd' w l,'X50.DB.BXZ.semd' w l
