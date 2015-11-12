
#set xrange [0:2000]

#plot "seismogram_F90.txt" w l lc 1, "seismogram_C_single_correct.txt" w l lc 3, "seismogram_C_single.txt" w l lc 4
plot "seismogram_C_single.txt" w l lc 1, "seismogram_C_single_correct.txt" w l lc 3

pause -1 "Hit any key..."

plot "seismogram_C_single.txt" w l lc 1
pause -1 "Hit any key..."

plot "seismogram_C_single_correct.txt" w l lc 3
pause -1 "Hit any key..."

