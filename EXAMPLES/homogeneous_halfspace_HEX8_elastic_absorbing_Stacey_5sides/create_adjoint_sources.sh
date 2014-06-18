#!/bin/bash

mkdir -p SEM

#cp OUTPUT_FILES/X20.DB.BX*.semd SEM/
if [ ! -e SEM/X20.DB.BXX.semd ]; then echo "please copy traces X20.DB.BX*.semd to SEM/"; exit 1; fi

#cp ~/SPECFEM3D_svn/utils/adjoint_sources/traveltime/xcreate_adjsrc_traveltime SEM/
if [ ! -e SEM/xcreate_adjsrc_traveltime ]; then echo "please make xcreate_adjsrc_traveltime and copy to SEM/"; exit 1; fi

# creates adjoint sources
cd SEM/

./xcreate_adjsrc_traveltime 10.0 25.0 3 X20.DB.BX*.semd
if [ ! -e X20.DB.BXZ.semd.adj ]; then echo "error creating adjoint sources, please check..."; exit 1; fi

rename .semd.adj .adj *semd.adj

cd ../

echo
echo

