#!/bin/csh

foreach file ($*)

set myfile = `basename $file .f90`

mv $file {$myfile}.f

end

