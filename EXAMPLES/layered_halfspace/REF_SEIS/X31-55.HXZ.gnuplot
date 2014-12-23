#!/usr/local/bin/gnuplot -persist
#
#    

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'DB.X31.HXZ.semd' w l,'DB.X55.HXZ.semd' w l
