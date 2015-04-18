
#set xrange [0:2000]

#plot "seismogram_F90.txt" w l 1, "seismogram_C_single_correct.txt" w l 3, "seismogram_C_single.txt" w l 4
plot "seismogram_C_single.txt" w l 1, "seismogram_C_single_correct.txt" w l 3

pause -1 "Hit any key..."

plot "seismogram_C_single.txt" w l 1
pause -1 "Hit any key..."

plot "seismogram_C_single_correct.txt" w l 3
pause -1 "Hit any key..."

