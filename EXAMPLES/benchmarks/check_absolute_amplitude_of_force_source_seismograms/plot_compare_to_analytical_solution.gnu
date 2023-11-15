set term x11

set xrange [0:1]

factor = 1

plot "Ux_time_analytical_solution_elastic.dat" w l lc 1, "CE.1.FXX.semd" using 1:(factor*$2) w l lc 3
pause -1 "Hit any key..."

plot "Uy_time_analytical_solution_elastic.dat" w l lc 1, "CE.1.FXY.semd" using 1:(factor*$2) w l lc 3
pause -1 "Hit any key..."

plot "Uz_time_analytical_solution_elastic.dat" w l lc 1, "CE.1.FXZ.semd" using 1:(factor*$2) w l lc 3
pause -1 "Hit any key..."

