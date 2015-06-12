#set size ratio -1
set term postscript eps font "Times-Roman,12"
set output "$1.eps"
#set noborder
#set noxtics ; set noytics
set title "receiver at "
plot "$1.dat" with lines
#, 'recfile_0005_disp_post_t.dat' with lines
#plot 'taup_prem_120.00km_rec7_0.deg_both.dat' using 1:1xticlabels(2) with linespoints
set xrange [ 0: $2]
set xlabel 'time [s]'
set ylabel 'displacement [m]' 
