#!/bin/bash

#############################
# First find hardpath files #
#############################

# Define rootdire
rootdir=`pwd`

# First lma's pathes then macro
MACRO_lma='/mnt/Data1/vmont/svn_lma/DSM_FOR_SPECFEM3D'
MACRO_macro='MY_ROOT_DIR_MACRO'

# Look for files to modify
find ./Examples/ -type f -exec grep -l $MACRO_lma {} \;   > my_founds_path_files_lma
find ./Examples/ -type f -exec grep -l $MACRO_macro {} \; > my_founds_path_files_macro

# Change lma ones
while read line
do
  sed -i "s|$MACRO_lma|$rootdir|g" $line
done < my_founds_path_files_lma

# Change macro ones
while read line
do
        sed -i "s|$MACRO_macro|$rootdir|g" $line
done < my_founds_path_files_macro


##### OLD
#MACRO='MY_ROOT_DIR_MACRO'
#find ./Examples/ -type f -exec sed -i "si|$MACRO|$rootdir|g" {} \;
# FIND HARDPATH Files
# find ./Examples/ -type f -exec grep -l /mnt/Data1 {} \;

