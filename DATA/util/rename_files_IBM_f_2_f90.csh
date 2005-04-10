#!/bin/csh

foreach file ($*)

set myfile = `basename $file .f`

mv $file {$myfile}.f90

end

