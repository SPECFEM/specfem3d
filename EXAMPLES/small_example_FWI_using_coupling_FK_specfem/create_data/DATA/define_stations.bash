# create grid of stations
# script to create station file for specfem
# compile tools in the folowing directory :
# (update compile.csh script, just fill the fortran compiler name)
dir=../../../../EXTERNAL_CODES_coupled_with_SPECFEM3D/AxiSEM_for_SPECFEM3D/utils_misc/


# numbers of stations
nb_sta_x=29
nb_sta_y=29

# distance betwee station (m)
distance_betxeen_sta_x=1000
distance_betxeen_sta_Y=1000

# coordinate of first station (m)
xmin=2000
ymin=2000

# station coordinate
depth=10000

#
NETWORK=XX
CODE=S

echo $nb_sta_x $nb_sta_y > Par_sta
echo $distance_betxeen_sta_x $distance_betxeen_sta_Y >> Par_sta
echo $xmin $ymin >> Par_sta
echo $depth >> Par_sta
echo $CODE >> Par_sta
echo $NETWORK >> Par_sta

$dir/create_uniform_stations.x<<EOF
Par_sta
EOF


cp station_list.txt STATIONS
