#!/bin/csh

set hdur = 11.2

foreach file ( $* )

set nlines = `wc -l $file `
echo $nlines > input_convolve_code.txt
echo $hdur >> input_convolve_code.txt
# use .true. for a triangle and .false. for a Gaussian
#echo ".true." >> input_convolve_code.txt
echo ".false." >> input_convolve_code.txt

echo convolving $file with hdur = $hdur using lines $nlines 

./xconvolve_source_timefunction < $file > ${file}.convolved
rm input_convolve_code.txt

end

