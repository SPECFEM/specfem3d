#!/usr/local/bin/gnuplot -persist
#
#    

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'X31.DB.HXZ.semd' w l,'X55.DB.HXZ.semd' w l
