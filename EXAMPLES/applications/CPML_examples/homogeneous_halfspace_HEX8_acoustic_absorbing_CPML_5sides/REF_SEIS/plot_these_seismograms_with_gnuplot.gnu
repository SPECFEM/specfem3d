#
# type " gnuplot plot_these_seismograms_with_gnuplot.gnu " to use this script
#    

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'DB.X20.MXZ.semd' w l,'DB.X30.MXZ.semd' w l,'DB.X40.MXZ.semd' w l,'DB.X50.MXZ.semd' w l
pause -1 "Hit key..."

