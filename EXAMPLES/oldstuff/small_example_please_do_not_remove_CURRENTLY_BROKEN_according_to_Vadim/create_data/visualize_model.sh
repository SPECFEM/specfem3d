#
#  USAGE
#
#   ./create_one_snaphot iteration
#
#   iteration is time step that you want to volume display
#
#

# specfem bin directory
bin=../../../specfem3d_Git_devel/bin

# choose what to display
DISPLAY=$1

# available choices :
#   velocity_X, velocity_Y, velocity_Z
#   curl_X, curl_Y, curl_Z
#   div_glob
#
#
#

# choose DATABASES_MPI DIRECTORY
DIRIN=DATABASES_MPI

# choose output directory
DIROUT=.

# choose resolution (low=0, high=1)
res=1

# --- DO NOT CHANGE------------------------------
declare -i NPROC
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NPROC="$NPROC-1"
it=$(printf "%06d" $1)
$bin/xcombine_vol_data_vtk 0 $NPROC $DISPLAY $DIRIN"/" $DIROUT $res
