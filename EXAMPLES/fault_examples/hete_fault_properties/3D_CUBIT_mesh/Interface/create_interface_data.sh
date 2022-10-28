###  Load Parameters file  ###
source ../Parameters_for_curved_surfaces.sh

###  Output
Subxyz=slabInterface_sigma_${Int_sigma}_inc_${Int_inc}.xyz

##-----------------------------------
##-----------------------------------
if [ ! -f ${Int_input_data} ] ; then
   echo "Miss the data for the curved fault. Please give the correct path to the data..."
   echo "Fails to create curved mesh..."
   exit
fi

if [ $Int_GRD_data = "1" ] ; then
# Read GRID-format data
  grd2xyz ${Int_input_data} | gawk '{if($1>= '"${X1}"' && $1<= '"${X2}"' && $2>= '"${Y1}"' && $2<='"${Y2}"') print ($1-('"$Lon_ref"'))*'"$Lon_scale"'*1e3, ($2-('"$Lat_ref"'))*'"$Lat_scale"'*1e3, -$3*1e3}' >  temp_file/output.xyz
else
# Read text data. Here is from Yue's slip model.
# For other slip model, one needs to change the number of head line (it's 9 lines here), and sort the order of columns in this sequence: longitude, latitude, depth.
  gawk '{if(NR>9 && $15==1) print $4,$3,$5}' ${Int_input_data} | blockmedian  -R${X1}/${X2}/${Y1}/${Y2} -I${Int_sample_inc} | sort -k1,1n -k2,2n > temp_file/output_mean.xyz
  # For xyz filei. Sort the order of columns in this sequence: longitude, latitude, depth.
  #blockmedian ${Int_input_data}  -R${X1}/${X2}/${Y1}/${Y2} -I${Int_sample_inc} | sort -k1,1n -k2,2n > temp_file/output_mean.xyz
  surface temp_file/output_mean.xyz  -Gtemp_file/output_surf.grd -R  -I${Int_sample_inc} -T1
  grd2xyz temp_file/output_surf.grd | gawk '{print ($1-('"$Lon_ref"'))*'"$Lon_scale"'*1e3, ($2-('"$Lat_ref"'))*'"$Lat_scale"'*1e3, -$3*1e3}' >  temp_file/output.xyz
fi

#####################################
#####      Smooth the data     ######
#####################################
$Run_python  ../python_scripts/Smooth_slab.py -p temp_file/output.xyz  ${Int_sigma} ${Int_inc}  > temp_file/${Subxyz}

#####################################
#####   Create jou file        ######
#####################################
dimension=`$Run_python ../python_scripts/Dimension.py -p  temp_file/${Subxyz}`
echo "Interface data dimension:" $dimension
$Run_python ../python_scripts/slab_interface_netsurf.py -p  temp_file/${Subxyz} ${dimension}   CUBIT_jou/slab_interface_netsurf_sigma_${Int_sigma}_inc_${Int_inc}.jou  ../output/interface_sigma_${Int_sigma}_inc_${Int_inc}.sat

#####################################
#####   Create sat file        ######
#####################################
$Run_python ../python_scripts/playback_jou.py -p  CUBIT_jou/slab_interface_netsurf_sigma_${Int_sigma}_inc_${Int_inc}.jou

#####################################
#####      Map to figure       ######
#####################################
final_Int_sample_inc=`gawk 'BEGIN{print '"$Int_sample_inc"'*'"${Int_inc}"'}'`
makecpt -Ctopo  -T0/50/5 -D -N > temp.cpt
gawk '{if(NR>9 && $15==1) print $4,$3,$5}' ${Int_input_data} | pscontour -R${X1}/${X2}/${Y1}/${Y2} -JX4i -B0.5f0.1:."Color Original; Lines Smoothed":WSne -I -Ctemp.cpt  -K > ps/${Int_sigma}_${Int_inc}_interface.eps
gawk '{print $1/1e3/'"$Lon_scale"'+('"$Lon_ref"'),$2/1e3/'"$Lat_scale"'+('"$Lat_ref"'),$3/-1e3}' temp_file/${Subxyz}  | xyz2grd -R -I${final_Int_sample_inc} -Gtemp_file/${Subxyz}.grd
grdcontour temp_file/${Subxyz}.grd -R -J -C5 -O  >> ps/${Int_sigma}_${Int_inc}_interface.eps


ps2pdf ps/${Int_sigma}_${Int_inc}_interface.eps ps/${Int_sigma}_${Int_inc}_interface.pdf
rm ps/${Int_sigma}_${Int_inc}_interface.eps temp.cpt  gmt.history -f

