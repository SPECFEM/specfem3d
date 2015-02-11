
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 ! number of processors =          168
 !
 ! number of ES nodes =    21.00000
 ! percentage of total 640 ES nodes =    3.281250      %
 ! total memory available on these ES nodes (Gb) =    336.0000
 !
 ! max points per processor = max vector length =       163941
 ! min vector length =           25
 ! min critical vector length =           75
 !
 ! on ES and SX-5, make sure "loopcnt=" parameter
 ! in Makefile is greater than       163941
 !
 ! total elements per AB slice =         2412
 ! total points per AB slice =       163941
 !
 ! total for full mesh:
 ! -------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !       405216
 ! approximate total number of points in entire mesh =
 !    27542088.0000000
 ! approximate total number of degrees of freedom in entire mesh =
 !    82626264.0000000
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along X =          336
 ! spectral elements along Y =          288
 ! GLL points along X =         1345
 ! GLL points along Y =         1153
 ! average distance between points along X in m =    475.4175
 ! average distance between points along Y in m =    436.8476
 !

 integer, parameter :: NSPEC_AB =         2412
 integer, parameter :: NGLOB_AB =       163941

 !
 ! number of time steps =        18200
 !
 integer, parameter :: NSPEC_ATTENUATION = NSPEC_AB
 logical, parameter :: ATTENUATION_VAL = .true.
 logical, parameter :: ANISOTROPY_VAL = .false.

