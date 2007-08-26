#!/bin/csh

# we mimic a triangle of half duration equal to half_duration_triangle
# using a Gaussian having a very close shape, as explained in Figure 4.2
# of the manual

set half_duration_triangle = 11.2

foreach file ( $* )

set nlines = `wc -l $file `
echo $nlines > input_convolve_code.txt
echo $half_duration_triangle >> input_convolve_code.txt
# use .true. for a triangle and .false. for a Gaussian
#echo ".true." >> input_convolve_code.txt
echo ".false." >> input_convolve_code.txt

echo convolving $file with half_duration_triangle = $half_duration_triangle using lines $nlines 

./xconvolve_source_timefunction < $file > ${file}.convolved

rm input_convolve_code.txt

end

