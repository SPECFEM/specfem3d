gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 6p LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 1 LABEL_FONT 1 HEADER_FONT_SIZE 12 FRAME_PEN 1p TICK_PEN 1p
makecpt -Cseis -T-1.010e-01/1.010e-01/3.108e-03 -D > color1a.cpt
makecpt -Cseis -T3.167e+00/3.868e+00/1.079e-02 -D > color1b.cpt
makecpt -Cseis -T-1.010e-01/1.010e-01/3.108e-03 -D > color2.cpt
psbasemap -JM2.8 -R-121.2/-114.8/32.3/36.7 -Ba1f0.5d:." ":wesn -K -V   -X1i -Y5i > horz_xc_vs_m16_m00_021_extra0.ps
psbasemap -JM2.8 -R-121.2/-114.8/32.3/36.7 -Ba1f0.5d:." ":wesn -G220/220/220 -K -O -V >> horz_xc_vs_m16_m00_021_extra0.ps
grep -v NaN INPUT/horz_02/horz_xc_vs_m16_m00_021_mask1.dat | awk '{print $1,$2,$3/1000}' | xyz2grd -Gtemp.grd -R-121.2/-114.8/32.3/36.7 -I2m
grdimage temp.grd -Ccolor1b.cpt -JM2.8 -Q -K -O -V >> horz_xc_vs_m16_m00_021_extra0.ps
psscale -Ccolor1b.cpt -D0.84/-0.2/1.68/0.15h -B2.00e-01f0.1:"  ": -E10p -K -O -V >> horz_xc_vs_m16_m00_021_extra0.ps
pstext -N -R0/1/0/1 -JM2.8 -K -O -V >>horz_xc_vs_m16_m00_021_extra0.ps<<EOF
0.7 -0.19 10 0 1 LM (3.50 \261 10 \045)
EOF
pstext -N -R0/1/0/1 -JM2.8 -K -O -V >>horz_xc_vs_m16_m00_021_extra0.ps<<EOF
0.7 -0.1 10 0 1 LM Vs  km/s
EOF
pscoast -JM2.8 -R-121.2/-114.8/32.3/36.7 -Df -A0/0/4 -W0.75p -Na/0.75p -K -O -V >> horz_xc_vs_m16_m00_021_extra0.ps
psxy /home/carltape/gmt/faults/jennings_more.xy -JM2.8 -R-121.2/-114.8/32.3/36.7 -m -W0.75p,0/0/0 -K -V -O >> horz_xc_vs_m16_m00_021_extra0.ps
psbasemap -JM2.8 -R-121.2/-114.8/32.3/36.7 -Ba1f0.5d:." ":wesn -K -V -O >> horz_xc_vs_m16_m00_021_extra0.ps
pstext -R0/1/0/1 -JM2.8 -N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0 -K -O -V >>horz_xc_vs_m16_m00_021_extra0.ps<<EOF
0 0.05 11 0 1 LB m00
EOF
pstext -N -R0/1/0/1 -JM2.8 -K -O -V >>horz_xc_vs_m16_m00_021_extra0.ps<<EOF
-0.15 0.4 14 90 1 CM Depth  =  10.0 km
EOF
psbasemap -JM2.8 -R-121.2/-114.8/32.3/36.7 -Ba1f0.5d:." ":WEsn -K -O -V -X3.3 >> horz_xc_vs_m16_m00_021_extra0.ps
psbasemap -JM2.8 -R-121.2/-114.8/32.3/36.7 -Ba1f0.5d:." ":WEsn -G220/220/220 -K -O -V >> horz_xc_vs_m16_m00_021_extra0.ps
grep -v NaN INPUT/horz_02/horz_xc_vs_m16_m00_021_mask1.dat | awk '{print $1,$2,$4/1000}' | xyz2grd -Gtemp.grd -R-121.2/-114.8/32.3/36.7 -I2m
grdimage temp.grd -Ccolor1b.cpt -JM2.8 -Q -K -O -V >> horz_xc_vs_m16_m00_021_extra0.ps
psscale -Ccolor1b.cpt -D0.84/-0.2/1.68/0.15h -B2.00e-01f0.1:"  ": -E10p -K -O -V >> horz_xc_vs_m16_m00_021_extra0.ps
pstext -N -R0/1/0/1 -JM2.8 -K -O -V >>horz_xc_vs_m16_m00_021_extra0.ps<<EOF
0.7 -0.19 10 0 1 LM (3.50 \261 10 \045)
EOF
pstext -N -R0/1/0/1 -JM2.8 -K -O -V >>horz_xc_vs_m16_m00_021_extra0.ps<<EOF
0.7 -0.1 10 0 1 LM Vs  km/s
EOF
pscoast -JM2.8 -R-121.2/-114.8/32.3/36.7 -Df -A0/0/4 -W0.75p -Na/0.75p -K -O -V >> horz_xc_vs_m16_m00_021_extra0.ps
psxy /home/carltape/gmt/faults/jennings_more.xy -JM2.8 -R-121.2/-114.8/32.3/36.7 -m -W0.75p,0/0/0 -K -V -O >> horz_xc_vs_m16_m00_021_extra0.ps
psbasemap -JM2.8 -R-121.2/-114.8/32.3/36.7 -Ba1f0.5d:." ":WEsn -K -V -O >> horz_xc_vs_m16_m00_021_extra0.ps
pstext -R0/1/0/1 -JM2.8 -N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0 -K -O -V >>horz_xc_vs_m16_m00_021_extra0.ps<<EOF
0 0.05 11 0 1 LB m16
EOF
psbasemap -JM2.8 -R-121.2/-114.8/32.3/36.7 -Ba1f0.5d:." ":wesn -K -O -V -X3.3 >> horz_xc_vs_m16_m00_021_extra0.ps
psbasemap -JM2.8 -R-121.2/-114.8/32.3/36.7 -Ba1f0.5d:." ":wesn -G220/220/220 -K -O -V >> horz_xc_vs_m16_m00_021_extra0.ps
awk '{print $1,$2,$5}' INPUT/horz_02/horz_xc_vs_m16_m00_021_mask1.dat | xyz2grd -Gtemp.grd -R-121.2/-114.8/32.3/36.7 -I2m
grdimage temp.grd -Ccolor2.cpt -JM2.8 -K -O -V -Q >> horz_xc_vs_m16_m00_021_extra0.ps
psscale -Ccolor2.cpt -D0.84/-0.2/1.68/0.15h -Ba1.00e-01f0.05:" ": -E10p -K -O -V >> horz_xc_vs_m16_m00_021_extra0.ps
pstext -N -R0/1/0/1 -JM2.8 -K -O -V >>horz_xc_vs_m16_m00_021_extra0.ps<<EOF
0.7 -0.1 10 0 1 LM ln(m16 / m00)
EOF
pscoast -JM2.8 -R-121.2/-114.8/32.3/36.7 -Df -A0/0/4 -W0.75p -Na/0.75p -K -O -V >> horz_xc_vs_m16_m00_021_extra0.ps
psxy /home/carltape/gmt/faults/jennings_more.xy -JM2.8 -R-121.2/-114.8/32.3/36.7 -m -W0.75p,0/0/0 -K -V -O >> horz_xc_vs_m16_m00_021_extra0.ps
psbasemap -JM2.8 -R-121.2/-114.8/32.3/36.7 -Ba1f0.5d:." ":wesn -K -V -O >> horz_xc_vs_m16_m00_021_extra0.ps
pstext -R0/1/0/1 -JM2.8 -N -C4p -W255/255/255o,1.0p,0/0/0,solid -G0/0/0 -K -O -V >>horz_xc_vs_m16_m00_021_extra0.ps<<EOF
0 0.05 11 0 1 LB ln(m16/m00)
EOF
pstext -N -R0/1/0/1 -JM2.8 -O -V >>horz_xc_vs_m16_m00_021_extra0.ps<<EOF
 0.5 0.98 16 0 1 CM
EOF
ps2pdf horz_xc_vs_m16_m00_021_extra0.ps
