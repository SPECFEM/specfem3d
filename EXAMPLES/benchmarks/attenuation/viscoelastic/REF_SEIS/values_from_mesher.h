
 !
 ! purely informative use
 !
 ! mesh statistics:
 ! ---------------
 !
 ! note:
 !    the values are only approximate and differ for different processes
 !    because the CUBIT + SCOTCH mesh has
 !    a different number of mesh elements and points in each slice
 !
 ! number of processors =           16
 !
 ! number of ES nodes =    2.000000
 ! percentage of total 640 ES nodes =   0.3125000      %
 ! total memory available on these ES nodes (Gb) =    32.00000
 !
 ! min vector length =           25
 ! min critical vector length =           75
 !
 ! master process: total points per AB slice =      2106081
 ! total elements per AB slice = (will be read in external file)
 ! total points per AB slice = (will be read in external file)
 !
 ! total for full mesh:
 ! -------------------
 !
 !
 ! number of time steps =         1300
 !
 ! time step =   5.000000000000000E-004
 !
 ! attenuation uses:
 !  NSPEC_ATTENUATION =        32000
 !
 ! anisotropy uses:
 !  NSPEC_ANISO =            1
 !
 ! adjoint uses:
 !  NSPEC_ADJOINT =            1
 !  NGLOB_ADJOINT =            1
 !
 ! approximate least memory needed by the solver:
 ! ----------------------------------------------
 !
 ! size of arrays for the largest slice =    218.360626220703       MB
 !                                      =   0.213242799043655       GB
 !
 !   (should be below 90% or so of the amount of memory available per processor
 core
 !   (if significantly more, the job will not run by lack of memory)
 !   (if significantly less, you waste a significant amount of memory)
 !
 ! check parameter to ensure the code has been compiled with the right values:
 &MESHER
 ABSORB_FREE_SURFACE_VAL = F
 /

