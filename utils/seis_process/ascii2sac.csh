#!/bin/csh

# checks arguments
if( $1 =~ "" ) then
  echo "Usage: ascii2sac.csh filenames"
  exit 1
endif

foreach file ($*)
  ./asc2sac $file
end
