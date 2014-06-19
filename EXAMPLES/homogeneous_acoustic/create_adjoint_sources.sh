#!/bin/bash

mkdir -p SEM

#cp OUTPUT_FILES/X20.DB.BX*.semd SEM/
if [ ! -e SEM/X20.DB.BXX.semd ]; then echo "please copy traces X20.DB.BX*.semd to SEM/"; exit 1; fi


#cp ~/SPECFEM3D_svn/utils/adjoint_sources/traveltime/xcreate_adjsrc_traveltime SEM/
if [ ! -e SEM/xcreate_adjsrc_traveltime ]; then echo "please make xcreate_adjsrc_traveltime and copy to SEM/"; exit 1; fi

# creates adjoint sources
cd SEM/

./xcreate_adjsrc_traveltime 9. 30. 3 X20.DB.BX*.semd
if [ ! -e X20.DB.BXZ.semd.adj ]; then echo "error creating adjoint sources, please check..."; exit 1; fi

rename .semd.adj .adj *semd.adj

# uses Z-component source for all 
# (acoustic adjoint sources are read from component 1 -> BXX.adj)
cp X20.DB.BXZ.adj X20.DB.BXX.adj
cp X20.DB.BXZ.adj X20.DB.BXY.adj

cd ../

echo
echo

