#!/bin/bash
#
#
# creates AVS_*.inp file for visualization (e.g. in Paraview)


./bin/xcreate_movie_shakemap_AVS_DX_GMT << EOF
2
1
2000
1
4
EOF

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo

