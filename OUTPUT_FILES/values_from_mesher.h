
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 ! number of processors =  152
 !
 ! number of ES nodes =  19.00000
 ! percentage of total 640 ES nodes =  2.968750  %
 ! total memory available on these ES nodes (Gb) =  304.0000
 !
 ! max points per processor = max vector length =  354777
 ! min vector length =  25
 ! min critical vector length =  75
 !
 ! on ES and SX-5, make sure "loopcnt=" parameter
 ! in Makefile is greater than  354777
 !
 ! total elements per AB slice =  5256
 ! total points per AB slice =  354777
 !
 ! total for full mesh:
 ! -------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !  798912
 ! approximate total number of points in entire mesh = 
 !  5.39261040000000E+7
 ! approximate total number of degrees of freedom in entire mesh = 
 !  1.61778312000000E+8
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along X =  456
 ! spectral elements along Y =  384
 ! GLL points along X =  1825
 ! GLL points along Y =  1537
 ! average distance between points along X in m =  314.0931
 ! average distance between points along Y in m =  329.3387
 !

 integer, parameter :: NSPEC_AB =  5256
 integer, parameter :: NGLOB_AB =  354777

 integer, parameter :: NSPEC_ATTENUATION = NSPEC_AB
 logical, parameter :: ATTENUATION_VAL = .true.
 logical, parameter :: ANISOTROPY_VAL = .false.

 integer, parameter :: SIMULATION_TYPE = 1
 logical, parameter :: SAVE_FORWARD = .false.


