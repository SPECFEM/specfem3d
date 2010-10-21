#!/bin/bash

gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 20 HEADER_FONT_SIZE 24

#=============== User input

# SPECFEM3D specific
dt=0.0005
t0=2.5
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

# Choose color scale
makecpt -Cjet -T0./1.5/.05 -V -Z > hor.cpt

# Loop on frame number (here the gmt filenames are gmt_movieH_??????.xyz)
for frame in `ls gmt_movieH*xyz |cut -c12-17`
do
      xyzfile="gmt_movieH_${frame}.xyz"
      grdfile=gmt_movieH_${frame}.grd
      psfile=snap_${frame}.ps
      jpgfile=snap_${frame}.jpg
      time=`echo "$dt * $frame -$t0" |bc -l `

psxy $JJJ $RRzoom -K -P -V -X4 -Y3 <<EOF > $psfile
EOF

xyz2grd $xyzfile -I$dx/$dy $RR -N0 -G$grdfile -V
grdimage $grdfile $JJ $RRzoom  -Chor.cpt -B0/0/10 -E50 -K -O -V -P >> $psfile
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
#$longmean $lattitle 32 0 1 CB Strong Motion Case 1
$longmean $lattitle 32 0 1 CB Horizontal Velocity
EOF

# Time step label
pstext $JJ $RRzoom -N -W255/255/255 -K -O -P -V <<EOF >>$psfile
$lonmaxzoom $latminzoom 24 0 1 RT Time = $time s
EOF

# Scale (to be set manually)
psscale -Chor.cpt -D${locscale}c/6c/8/.5 -O -B${scaleincrement}::/:${scaleunit}: <<EOF >>$psfile
EOF

rm -f $grdfile
convert -quality 100 $psfile $jpgfile
#display $psfile
rm -f $psfile
done

