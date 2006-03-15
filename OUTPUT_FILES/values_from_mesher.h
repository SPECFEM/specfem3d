 
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 ! number of processors =           36
 !
 ! number of ES nodes =    4.500000    
 ! percentage of total 640 ES nodes =   0.7031250      %
 ! total memory available on these ES nodes (Gb) =    72.00000    
 !
 ! max points per processor = max vector length =       703989
 ! min vector length =           25
 ! min critical vector length =           75
 !
 ! on ES and SX-5, make sure "loopcnt=" parameter
 ! in Makefile is greater than       703989
 !
 ! total elements per AB slice =        10512
 ! total points per AB slice =       703989
 !
 ! total for full mesh:
 ! -------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !       378432
 ! approximate total number of points in entire mesh = 
 !    25343604.0000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    76030812.0000000     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along X =          288
 ! spectral elements along Y =          288
 ! GLL points along X =         1153
 ! GLL points along Y =         1153
 ! average distance between points along X in m =    239.0178    
 ! average distance between points along Y in m =    287.6366    
 !
 
 integer, parameter :: NSPEC_AB =        10512
 integer, parameter :: NGLOB_AB =       703989
 
 integer, parameter :: NSPEC_ATTENUATION = 1
 logical, parameter :: ATTENUATION_VAL = .false.
 logical, parameter :: ANISOTROPY_VAL = .false.
 
