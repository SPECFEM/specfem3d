#!/bin/csh

foreach file ($*)
  ./asc2sac $file
end
