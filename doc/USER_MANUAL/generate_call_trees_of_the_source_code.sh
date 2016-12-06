#!/bin/sh

# changes to wiki/ directory
if [ -d ../call_trees_of_the_source_code ]; then
  cd ../call_trees_of_the_source_code/
else
  echo "call_trees_of_the_source_code/ directory does not exist, nothing to do..."
  exit 1
fi

# runs script
./run_doxygen.sh

echo
echo

