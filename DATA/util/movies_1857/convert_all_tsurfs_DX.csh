#!/bin/csh

# convert all GoCaD t_surfs to OpenDX

# call with " convert_all_tsurf_DX.csh *.ts "

rm input.ts

foreach file ( $* )

echo converting $file to OpenDX
echo number of TFACE in file is `grep TFACE $file | wc -l`

cp $file input.ts
./xdraw
set newfile = `basename $file .ts`
mv DX_fullmesh.dx DX_{$newfile}.dx
# DK DK also save in ASCII format
mv ASCII_fullmesh.dat ASCII_{$newfile}.dat

echo " "

end

rm input.ts

