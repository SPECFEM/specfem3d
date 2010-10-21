#!/bin/bash

gmtset BASEMAP_TYPE fancy ANOT_FONT_SIZE 20 HEADER_FONT_SIZE 24 BASEMAP_AXES WeSn

#=============== User input

# SPECFEM3D specific
lonmin=5.6
lonmax=6.0
latmin=44.98
latmax=45.25
nlon=193
nlat=193

# Plot specific
lonminzoom=5.6
lonmaxzoom=6.0
latminzoom=45.05
latmaxzoom=45.25
sizeplot=25
scaleunit="m/s"
scaleincrement="0.2"
# Choose color scale
makecpt -Cjet -T0./1./.05 -V -Z > pgv.cpt

# Give PGV filenames
   xyzfile="gmt_shaking_map_Grenoble_basin_vel.geo"
   grdfile="gmt_shaking_map_Grenoble_basin_vel.grd"
   psfile="gmt_shaking_map_Grenoble_basin_vel.ps"
   jpgfile="gmt_shaking_map_Grenoble_basin_vel.jpg"

# Grenoble basin specific
case="S1"

#=============== End of user input

dx=`echo "($lonmax - $lonmin ) / ($nlon -1)" |bc -l `
dy=`echo "($latmax - $latmin )/($nlat -1)" |bc -l `
bigsize=`echo "$sizeplot * 1.15 " |bc -l `
locscale=`echo "$sizeplot * 1.05 " |bc -l `
longmean=`echo "($lonminzoom + $lonmaxzoom) / 2 " |bc -l `
lattitle=`echo "$latmaxzoom + .01 " |bc -l `
echo "title set at $longmean $lattitle "
JJ="-JM${sizeplot}c"
JJJ="-JM${bigsize}c"
RR="-R${lonmin}/${lonmax}/${latmin}/${latmax}"
RRzoom="-R${lonminzoom}/${lonmaxzoom}/${latminzoom}/${latmaxzoom}"

psxy $JJJ $RRzoom -K -P -V -X4 -Y3 <<EOF > $psfile
EOF

xyz2grd $xyzfile -I$dx/$dy $RR -N0 -G$grdfile -V
grdimage $grdfile $JJ $RRzoom  -Cpgv.cpt -B5m/5m/10 -E50 -K -O -V -P >> $psfile
# Contourfile "contour_hor.dat" has to be defined first
grdcontour -S -Q2 $grdfile $JJ $RRzoom -W3/0/0/0 -Ccontour_hor.dat -K -O -V -P >> $psfile
# Grenoble basin specific
grdcontour bedrock_map_GEO.grd $JJ $RRzoom -W5/255/255/255 -Ccontour_bassin.dat -K -O -V -P >>$psfile

#------------- Grenoble basin specific labels
#Source
if [ ${case} == "W1" ]
then
psxy $JJ $RRzoom -Sa.5 -W10 -K -O -P -V  <<EOF >> $psfile
5.9167 45.2167
EOF
fi
if [ ${case} == "W2" ]
then
psxy $JJ $RRzoom -Sa.5 -W10 -K -O -P -V  <<EOF >> $psfile
5.75 45.025
EOF
fi
if [ ${case} == "S1" ]
then
psxy $JJ $RRzoom -W12/255/255/255 -K -O -P -V  <<EOF >> $psfile
5.875621 45.189701
5.957818 45.243744
EOF
fi
if [ ${case} == "S2" ]
then
psxy $JJ $RRzoom -W12/255/255/255 -K -O -P -V  <<EOF >> $psfile
5.797558 45.003975
5.702406 45.046051
EOF
fi
#------------- End of Grenoble basin specific labels

# Title of plot
pstext $JJ $RRzoom -N -W255/255/255 -K -O -P -V -K <<EOF >>$psfile
$longmean $lattitle 32 0 1 CB PGV (m/s)
EOF

# Scale (to be set manually)
psscale -Cpgv.cpt -D${locscale}c/6c/8/.5 -O -B${scaleincrement}::/:${scaleunit}: <<EOF >>$psfile
EOF

rm -f $grdfile

# use ImageMagick to convert the file to JPEG with maximum quality
convert -quality 100 $psfile $jpgfile

#display $psfile

rm -f $psfile

