#
# type " gnuplot plot_these_seismograms_with_gnuplot.gnu " to use this script
#    

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'X20.DB.BXZ.semd' w l,'X30.DB.BXZ.semd' w l,'X40.DB.BXZ.semd' w l,'X50.DB.BXZ.semd' w l
pause -1 "Hit key..."

