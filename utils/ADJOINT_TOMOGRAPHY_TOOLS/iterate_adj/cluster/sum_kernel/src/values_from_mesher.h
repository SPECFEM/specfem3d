
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
 integer, parameter ::  NPROC_VAL =          168
 integer, parameter :: NPROC_XI_VAL =           14
 integer, parameter :: NPROC_ETA_VAL =           12
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

 integer, parameter :: NSPEC_AB_VAL =         2412
 integer, parameter :: NGLOB_AB_VAL =       163941

 !
 ! number of time steps =        18200
 !
 integer, parameter :: NSPEC_ATTENUATION =            1
 logical, parameter :: ATTENUATION_VAL = .false.

 integer, parameter :: NSPEC_ANISO =            1
 logical, parameter :: ANISOTROPY_VAL = .false.

 integer, parameter :: NSPEC_ATT_AND_KERNEL =            1
 integer, parameter :: NSPEC_ADJOINT =         2412
 integer, parameter :: NGLOB_ADJOINT =       163941

 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_VAL =          150
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_VAL =          150
 integer, parameter :: NSPEC2D_BOTTOM_VAL =           36
 integer, parameter :: NSPEC2D_TOP_VAL =          576
 integer, parameter :: NPOIN2DMAX_XMIN_XMAX_VAL =         3751
 integer, parameter :: NPOIN2DMAX_YMIN_YMAX_VAL =         3751
 integer, parameter :: NPOIN2DMAX_XY_VAL =         3751

 integer, parameter :: NSPEC2D_MOHO_BOUN =            1
 integer, parameter :: NSPEC_BOUN =            1
