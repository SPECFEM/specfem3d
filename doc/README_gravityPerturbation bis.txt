SPECFEM3D has the capability to compute the perturbations of gravity 
induced by seismic waves at arbitrary locations above the Earth's surface. 
The current file documents the usage of this feature and the related code 
modifications.


Theory
------

The computation of gravity perturbations induced by deformation is based on 
equation 4 of:
J. Harms, R. DeSalvo, S. Dorsher and V. Mandic (2009), "Simulation of underground
gravity gradients from stochastic seismic fields", Phys. Rev. D, 80, 122001


Enabling gravity computations
-----------------------------

Gravity perturbation can be computed in any SPECFEM3D simulation by 
placing a file called "gravity_stations" in the DATA directory, in 
addition to the regular input files.


Input file format
-----------------

The format of the "gravity_stations" input file is:

n	dt_gap
x1	y1	z1
x2	y2	z2
...	...	...
xn	yn 	zn

where
 n : number of stations where gravity time series are needed
 dt_gap : gravity time series are sampled every dt_gap time steps of the 
          SPECFEM3D simulation


Output file format
-----------------

Time series of gravity are output in files named "OUTPUT_FILES/stat*.grav",
where * is the station index (one file per station). Their format is four
columns:
  t  ax  ay  az 
(time and acceleration along x, y and z, respectively).


Code modifications
------------------

All the routines related to the gravity perturbation computations are 
placed in one module called specfem3d/gravity_perturbation.f90. The module 
provides three public subroutines and a flag. Each subroutine is invoked 
in one of the following stages:
 1. during the initialization, 
 2. during the iterative time stepping scheme and 
 3. at the output stage.
In #1 the code checks for the presence of the input file "gravity_stations". 
If the file exists, the flag "GRAVITY_SIMULATION" is turned on and the 
subroutines #2 and #3 are invoked.


Author
------

Surendra Somala (Caltech) surendra@caltech.edu - 2013
with advice from Jan Harms and Pablo Ampuero (ampuero@gps.caltech.edu)
