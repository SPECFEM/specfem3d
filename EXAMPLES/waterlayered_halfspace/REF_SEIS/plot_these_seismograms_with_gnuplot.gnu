#
# type " gnuplot plot_these_seismograms_with_gnuplot.gnu " to use this script
#    

set xlabel "time (s)"
set ylabel "displacement (m)"

plot 'DB.X10.HXZ.semd' w l,'DB.X20.HXZ.semd' w l,'DB.X30.HXZ.semd' w l,'DB.X40.HXZ.semd' w l,'DB.X50.HXZ.semd' w l,'DB.X60.HXZ.semd' w l,'DB.X70.HXZ.semd' w l,'DB.X80.HXZ.semd' w l,'DB.X90.HXZ.semd' w l

pause -1 "Hit key..."

