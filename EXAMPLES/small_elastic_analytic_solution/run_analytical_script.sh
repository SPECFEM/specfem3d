#!/bin/bash

echo
echo "creating analytic solutions..."
echo

./analytical_solution_Aki.py 3000.0 1000.0 4000.0 DB.Z1.FX
./analytical_solution_Aki.py 1000.0 3000.0 4000.0 DB.Z2.FX
./analytical_solution_Aki.py 6000.0 1000.0 4000.0 DB.Z3.FX
./analytical_solution_Aki.py 1000.0 6000.0 4000.0 DB.Z4.FX
./analytical_solution_Aki.py 6000.0 1000.0 2000.0 DB.Z5.FX
./analytical_solution_Aki.py 1000.0 6000.0 2000.0 DB.Z6.FX

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mv -v  DB.Z* REF_ANALYTIC/

echo
echo "done"
echo
