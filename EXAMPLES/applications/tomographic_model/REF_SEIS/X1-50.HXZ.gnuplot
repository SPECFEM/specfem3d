#!/usr/bin/gnuplot -persist

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'DB.X1.HXZ.semd' w l,'DB.X10.HXZ.semd' w l,'DB.X20.HXZ.semd' w l,'DB.X30.HXZ.semd' w l,'DB.X40.HXZ.semd' w l,'DB.X50.HXZ.semd' w l

