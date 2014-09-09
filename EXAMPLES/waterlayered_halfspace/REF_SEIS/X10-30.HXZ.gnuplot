#!/data1/dpeter/install//bin/gnuplot -persist

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'X10.DB.HXZ.semd' w l,'X20.DB.HXZ.semd' w l,'X30.DB.HXZ.semd' w l

