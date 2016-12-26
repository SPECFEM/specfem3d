To generate the subduction zone mesh for SPECFEM3D:

1. In the linux shell, run:
	python process_slab_rotate.py
   This step creates a file called slab_rotate_before_loft.cub
2. Open the Trelis graphic window, then run create_chunks_mesh.jou
   This step creates a file called slab_rotate.cub
(2).
If you have no easy access to Trelis graphics window, you can run “trelis -nographics -input ./create_chunks_mesh.jou” instead.
 
3. In the linux shell, run: 
	python exportmesh.py 
   This step creates the SPECFEM3D mesh files

Notes:
Before starting, include the path to cubit.py in your $PYTHONPATH environment variable. 
Python scripts have been tested with python 2.7, but 2.6 to 3.5 should also be fine. 

Step 1 runs on the linux shell, but not directly in the Trelis window (the numpy and matplotlib modules cannot be found).
type "sh run.sh" to run all three steps.
