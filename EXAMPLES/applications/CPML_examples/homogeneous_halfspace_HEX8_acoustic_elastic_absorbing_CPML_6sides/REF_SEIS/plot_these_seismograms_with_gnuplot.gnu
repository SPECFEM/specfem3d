#
# type " gnuplot plot_these_seismograms_with_gnuplot.gnu " to use this script
#    

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'DB.X1.BXZ.semd' w l,'DB.X10.BXZ.semd' w l,'DB.X20.BXZ.semd' w l,'DB.X30.BXZ.semd' w l,'DB.X40.BXZ.semd' w l,'DB.X50.BXZ.semd' w l
pause -1 "Hit key..."

