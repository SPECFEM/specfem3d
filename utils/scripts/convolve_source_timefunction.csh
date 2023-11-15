#!/bin/csh

# We mimic a triangle of half duration equal to half_duration_triangle
# using a Gaussian having a very close shape, as explained in Figure 5.2
# of the manual

set half_duration_triangle = 1  # TODO : It shall be better to send this to the script with an option. Ex : ./convolve_source_timefunction fileToConvolve --hdur 11.2
set SPECFEM_PATH = "~/specfem3d"  # Path to your specfem3d/ directory. TODO Could also be send as an (optional) option --path_to_specfem
set use_triangle_source = ".false." # Use .true. for a triangle and .false. for a Gaussian. TODO Idem --use_triangle_source

########### DO NOT CHANGE ANYTHING BELOW ###########

foreach file ( $* ) # Loop on all files (and directories) given in argument

  set nlines = `wc -l $file`                              # nlines is the number of lines in the file considered
  echo $nlines > input_convolve_code.txt                  # Create a new temp file input_convolve_code.txt (or recreate it if it exists) and write nlines on it
  echo $half_duration_triangle >> input_convolve_code.txt # Write half_duration_triangle after that...
  echo $use_triangle_source >> input_convolve_code.txt    #  ... and write .true. or false at the end depending on the source chosen
  echo convolving $file with half_duration_triangle = $half_duration_triangle using lines $nlines # A small print
  $SPECFEM_PATH/bin/xconvolve_source_timefunction < $file > ${file}.convolved                     # Execute xconvolve_source_timefunction with argument $file and write the output on ${file}.convolved
  rm input_convolve_code.txt                              # Remove the temp file input_convolve_code.txt

end # End of the loop
