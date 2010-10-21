#!/usr/local/bin/gnuplot -persist
#
#    

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'X1.DB.BHZ.semd' w l,'X10.DB.BHZ.semd' w l,'X20.DB.BHZ.semd' w l,'X30.DB.BHZ.semd' w l,'X40.DB.BHZ.semd' w l,'X50.DB.BHZ.semd' w l
