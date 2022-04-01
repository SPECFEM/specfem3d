###  Load Parameters file  ###
source ../Parameters_for_curved_surfaces.sh

##################################
##################################
###  Output
Subxyz=surface_sigma_${Sur_sigma}_inc_${Sur_inc}.xyz


if [ $Sur_GRD_data = "1" ] ; then
    if [ ! -f $Sur_input_data ] ; then
       echo "Miss the data for the curved free surface. Please give the correct path to the data..."
       echo "Fails to create curved mesh..."
       exit
    fi
    # Because the global topography data may be very large, it is better to cut the global data to smaller one.
    if [ ! -f temp_file/cutted_${X1}_${X2}_${Y1}_${Y2}.grd ] ; then
       grdcut $Sur_input_data -Gtemp_file/cutted_${X1}_${X2}_${Y1}_${Y2}.grd -R${X1}/${X2}/${Y1}/${Y2}
    fi
    grd2xyz temp_file/cutted_${X1}_${X2}_${Y1}_${Y2}.grd | gawk '{print ($1-('"$Lon_ref"'))*'"$Lon_scale"'*1e3, ($2-('"$Lat_ref"'))*'"$Lat_scale"'*1e3, $3}' >  temp_file/topo.xyz
else
    gawk '{if($1>= '"${X1}"' && $1<= '"${X2}"' && $2>= '"${Y1}"' && $2<='"${Y2}"') print ($1-('"$Lon_ref"'))*'"$Lon_scale"'*1e3, ($2-('"$Lat_ref"'))*'"$Lat_scale"'*1e3, $3}' ${input_data} > temp_file/topo.xyz
fi

#####################################
#####      Smooth the data     ######
#####################################
$Run_python ../python_scripts/Smooth_slab.py -p  temp_file/topo.xyz  $Sur_sigma $Sur_inc  > temp_file/${Subxyz}

#####################################
#####   Create jou file        ######
#####################################
dimension=`$Run_python  ../python_scripts/Dimension.py -p  temp_file/${Subxyz}`
echo "Topography data dimension:" $dimension
$Run_python ../python_scripts/slab_interface_netsurf.py -p  temp_file/${Subxyz} ${dimension}   CUBIT_jou/surface_sigma_${Sur_sigma}_inc_${Sur_inc}.jou  ../output/surface_sigma_${Sur_sigma}_inc_${Sur_inc}.sat

#####################################
#####   Create sat file        ######
#####################################
$Run_python ../python_scripts/playback_jou.py -p  CUBIT_jou/surface_sigma_${Sur_sigma}_inc_${Sur_inc}.jou

#####################################
#####      Map to figure       ######
#####################################
makecpt -Crelief  -T-8000/8000/100 -N > temp.cpt
### Original data
grdimage   temp_file/cutted_${X1}_${X2}_${Y1}_${Y2}.grd   -R${X1}/${X2}/${Y1}/${Y2} -JX4i -B0.5f0.1/0.5f0.1:."Original topography":WSne -Ctemp.cpt -K -Q > ps/${Sur_sigma}_${Sur_inc}_topography.eps
grdcontour temp_file/cutted_${X1}_${X2}_${Y1}_${Y2}.grd   -R -J -C500 -A- -O -K >> ps/${Sur_sigma}_${Sur_inc}_topography.eps
### Smoothed data
LowerX=`gawk 'BEGIN{print ('"$X1"'-'"$Lon_ref"')*'"$Lon_scale"'*1e3}'`
UpperX=`gawk 'BEGIN{print ('"$X2"'-'"$Lon_ref"')*'"$Lon_scale"'*1e3}'`
LowerY=`gawk 'BEGIN{print ('"$Y1"'-'"$Lat_ref"')*'"$Lat_scale"'*1e3}'`
UpperY=`gawk 'BEGIN{print ('"$Y2"'-'"$Lat_ref"')*'"$Lat_scale"'*1e3}'`
pscontour temp_file/${Subxyz}  -R${LowerX}/${UpperX}/${LowerY}/${UpperY} -B50000f10000:."Smoothed and resampled topography":WSne -J -Ctemp.cpt -I -A- -O -K -X5i  >> ps/${Sur_sigma}_${Sur_inc}_topography.eps
pscontour temp_file/${Subxyz}  -R${LowerX}/${UpperX}/${LowerY}/${UpperY}  -J -C500 -Wthinnest -A- -O   >> ps/${Sur_sigma}_${Sur_inc}_topography.eps


ps2pdf ps/${Sur_sigma}_${Sur_inc}_topography.eps ps/${Sur_sigma}_${Sur_inc}_topography.pdf
rm ps/${Sur_sigma}_${Sur_inc}_topography.eps temp.cpt  gmt.history -f


