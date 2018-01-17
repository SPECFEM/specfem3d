#!/bin/bash

currentdir=`pwd`

# clean
rm -rf html/

# generates html files
doxygen Doxyfile


echo
echo "done"
echo
echo "for documentation, see:"
echo "> open $currentdir/html/index.html"
echo

