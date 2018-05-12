set term x11

set xrange [0:1]

plot "Ux_time_analytical_solution_elastic.dat" w l lc 1, "Ux_time_analytical_solution_viscoelastic.dat" w l lc 3, "Ux_time_analytical_solution_viscoelastic_without_near_field.dat" w l lc 4
pause -1 "Hit any key..."

plot "Uy_time_analytical_solution_elastic.dat" w l lc 1, "Uy_time_analytical_solution_viscoelastic.dat" w l lc 3, "Uy_time_analytical_solution_viscoelastic_without_near_field.dat" w l lc 4
pause -1 "Hit any key..."

plot "Uz_time_analytical_solution_elastic.dat" w l lc 1, "Uz_time_analytical_solution_viscoelastic.dat" w l lc 3, "Uz_time_analytical_solution_viscoelastic_without_near_field.dat" w l lc 4
pause -1 "Hit any key..."

