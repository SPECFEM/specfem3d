###  Setup reference point and range ####
# This parameters are only for the curved free surface or fault
# X1 and X2 are the lower and upper range of longitude
# Y1 and Y2 are the lower and upper range of latitude
# Lon_ref and Lat_ref are the reference point (i.e., (0,0) in Cartesian coordinates)
# It is a good choice to set up it as the epicenter.
X1=129.0
X2=133.5
Y1=31.2
Y2=34.7
Lon_ref=130.76
Lat_ref=32.76

# Python path
# Or set up the full path for python manually in your system
Run_python=`which python`
Work_dir="/u/moana/user/weng/Weng/Scripts_dev/3D_CUBIT_mesh"

###  Input data for curved surface
#GRD_data=0           # do not use grid data
Sur_GRD_data=1           # do use grid data
Sur_input_data="${Work_dir}/Surface/data/Local_topo.grd"
###  Smooth parameters
# Resample interval
Sur_inc=12
# Smooth parameter
Sur_sigma=1

###  Input data for curved surface
Int_GRD_data=0           # do not use grid data
#GRD_data=1           # do use grid data
Int_input_data="${Work_dir}/Interface/data/Slip_model_kumamoto.dat"
###  Smooth parameters
Int_inc=12
Int_sigma=1
Int_sample_inc=0.007



### Don't need to change
Earth_radius="6371.0"
Lon_scale=`gawk 'BEGIN{print cos('"$Lat_ref"'/180.0*3.1415926)*3.1415926*'"$Earth_radius"'/180.0}'`
Lat_scale=`gawk 'BEGIN{print 3.1415926*'"$Earth_radius"'/180.0}'`

