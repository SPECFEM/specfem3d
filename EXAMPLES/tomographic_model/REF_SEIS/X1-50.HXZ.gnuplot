#!/usr/bin/gnuplot -persist

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'X1.DB.HXZ.semd' w l,'X10.DB.HXZ.semd' w l,'X20.DB.HXZ.semd' w l,'X30.DB.HXZ.semd' w l,'X40.DB.HXZ.semd' w l,'X50.DB.HXZ.semd' w l

