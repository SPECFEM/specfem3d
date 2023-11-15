#!/bin/csh -f

# DK DK filter and convert Northridge data to velocity in SAC
# DK DK read sac ALPHA (ascii) format as input

rm -r xconv

ifort -o xconv convert_files_sacalpha_ascii.f90

foreach file ($*)
  echo $file
  set newfile = `basename $file .dataa_ascii.sac.veloc`
  ./xconv < ${file} > ${newfile}.datav
end

