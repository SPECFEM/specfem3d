#!/data1/dpeter/install//bin/gnuplot -persist

set xlabel "time (s)" 
set ylabel "displacement (m)" 

plot 'K400.CE.BXZ.semd' w l,'K851.CE.BXZ.semd' w l,'K853.CE.BXZ.semd' w l

