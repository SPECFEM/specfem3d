To generate a subduction zone mesh with Trelis for dynamic rupture simulations in SPECFEM3D:

0. Include the path to cubit.py in your $PYTHONPATH environment variable.
   Create in this directory a file "etopo2.xyz" containing the topography of the region of interest
   or a global topography file, e.g. https://www.ngdc.noaa.gov/mgg/global/etopo2.html
   (3 columns: longitude, latitude, elevation in meters)
   Create in this directory a file containing the fault geometry
   (3 columns: longitude, latitude, elevation in km)
   Modify the user parameters in file process_slab_rotate.py and exportmesh.py following the guidelines therein

1. In the linux shell, run:
     python process_slab_rotate.py
   This step creates a file called slab_rotate_before_loft.cub

2. Open the Trelis graphical user interface, then run create_chunks_mesh.jou
   This step creates a file called slab_rotate.cub
   If you have no easy access to the Trelis GUI, you can instead run in the linux shell:
     trelis -nographics -input ./create_chunks_mesh.jou

3. In the linux shell, run:
     python exportmesh.py
   This step creates the SPECFEM3D mesh files

Notes:
Python scripts have been tested with python 2.7, but 2.6 to 3.5 should also be fine.
Step 1 runs on the linux shell, but not directly in the Trelis window (the numpy and matplotlib modules cannot be found).
Type "sh run.sh" to run steps 1 to 3 automatically.
