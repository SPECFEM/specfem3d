#!/bin/csh -f
#
# Vala Hjorleifsdottir and Qinya Liu, Caltech, Jan 2007
#
# uses `asc2sac` binary which can be compiled from UTILS/lib/asc2sac.c
# (requires SAC libraries for compilation)

foreach file ($*)
	echo $file
  set nlines = `wc -l $file | awk '{print $1}'`
  /opt/seismo-util/bin/asc2sac $file $nlines $file.sac
end
