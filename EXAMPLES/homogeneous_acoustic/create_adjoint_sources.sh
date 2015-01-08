#!/bin/bash

mkdir -p SEM

# station
station="X20"
network="DB"

sta=$network.$station

#cp OUTPUT_FILES/$sta.BX*.semd SEM/
if [ ! -e SEM/$sta.BXX.semd ]; then echo "please copy traces $sta.BX*.semd to SEM/"; exit 1; fi


#cp ~/SPECFEM3D_svn/utils/adjoint_sources/traveltime/xcreate_adjsrc_traveltime SEM/
if [ ! -e SEM/xcreate_adjsrc_traveltime ]; then echo "please make xcreate_adjsrc_traveltime and copy to SEM/"; exit 1; fi

# creates adjoint sources
cd SEM/

./xcreate_adjsrc_traveltime 9. 30. 3 $sta.BX*.semd
if [ ! -e $sta.BXZ.semd.adj ]; then echo "error creating adjoint sources, please check..."; exit 1; fi

rename .semd.adj .adj *semd.adj

# uses Z-component source for all 
# (acoustic adjoint sources are read from component 1 -> BXX.adj)
cp $sta.BXZ.adj $sta.BXX.adj
cp $sta.BXZ.adj $sta.BXY.adj

cd ../

echo
echo

