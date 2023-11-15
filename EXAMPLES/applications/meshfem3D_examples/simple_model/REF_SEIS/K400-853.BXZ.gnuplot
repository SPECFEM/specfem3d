#!/data1/dpeter/install//bin/gnuplot -persist

set xlabel "time (s)" 
set ylabel "displacement (m)" 

plot 'CE.K400.BXZ.semd' w l,'CE.K851.BXZ.semd' w l,'CE.K853.BXZ.semd' w l

