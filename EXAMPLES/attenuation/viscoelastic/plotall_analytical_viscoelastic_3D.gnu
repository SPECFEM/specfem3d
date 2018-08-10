set term x11

set xrange [0.1:0.5]

plot "Ux_time_analytical_solution_elastic.dat" w l lc 1, "Ux_time_analytical_solution_viscoelastic.dat" w l lc 3, "CE.1.FXX.semd" title 'SPECFEM Ux with viscoelasticity' w l lc 5
pause -1 "Hit any key..."

plot "Uy_time_analytical_solution_elastic.dat" w l lc 1, "Uy_time_analytical_solution_viscoelastic.dat" w l lc 3, "CE.1.FXY.semd" title 'SPECFEM Uy with viscoelasticity' w l lc 5
pause -1 "Hit any key..."

plot "Uz_time_analytical_solution_elastic.dat" w l lc 1, "Uz_time_analytical_solution_viscoelastic.dat" w l lc 3, "CE.1.FXZ.semd" title 'SPECFEM Uz with viscoelasticity' w l lc 5
pause -1 "Hit any key..."
