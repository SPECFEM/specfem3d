#!/bin/csh

foreach file ($*)
	echo $file
  set nlines = `wc -l $file | awk '{print $1}'`
  ./asc2sac/asc2sac $file $nlines $file.sac
end
