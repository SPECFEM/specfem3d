
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 !
 ! number of chunks =            1
 !
 ! these statistics do not include the central cube
 !
 ! number of processors =          100
 !
 ! maximum number of points per region =       475797
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      1427391
 !
 ! total elements per slice =         7821
 ! total points per slice =       528135
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh =
 !    4692600.00000000
 ! approximate total number of points in entire mesh =
 !    316881000.000000
 ! approximate total number of degrees of freedom in entire mesh =
 !    892096200.000000
 !
 ! position of the mesh chunk at the surface:
 ! -----------------------------------------
 !
 ! angular size in first direction in degrees =    68.00000
 ! angular size in second direction in degrees =    68.00000
 !
 ! longitude of center in degrees =    7.000000
 ! latitude of center in degrees =    52.00000
 !
 ! angle of rotation of the first chunk =    45.00000
 !
 ! corner            1
 ! longitude in degrees =    6.99999999999999
 ! latitude in degrees =    8.40726023480213
 !
 ! corner            2
 ! longitude in degrees =    64.1611658328499
 ! latitude in degrees =    34.9449365466497
 !
 ! corner            3
 ! longitude in degrees =   -50.1611658328499
 ! latitude in degrees =    34.9449365466497
 !
 ! corner            4
 ! longitude in degrees =   -173.000000000000
 ! latitude in degrees =    84.3892903527184
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =          960
 ! GLL points along a great circle =         3840
 ! average distance between points in degrees =   9.3750000E-02
 ! average distance between points in km =    10.42452
 ! average size of a spectral element in km =    41.69810
 !
 ! number of time steps =        10600
 !
 ! number of seismic sources =            1
 !

 ! approximate static memory needed by the solver:
 ! ----------------------------------------------
 !
 ! size of static arrays per slice =   0.260327711701393       GB
 !
 !   (should be below and typically equal to 80% or 90%
 !    of the memory installed per core)
 !   (if significantly more, the job will not run by lack of memory)
 !   (if significantly less, you waste a significant amount of memory)
 !
 ! size of static arrays for all slices =    26.0327711701393       GB
 !                                      =   2.542262809583917E-002  TB
 !

 integer, parameter :: NEX_XI_VAL =          240
 integer, parameter :: NEX_ETA_VAL =          240

 integer, parameter :: NSPEC_CRUST_MANTLE =         7092
 integer, parameter :: NSPEC_OUTER_CORE =          684
 integer, parameter :: NSPEC_INNER_CORE =           45

 integer, parameter :: NGLOB_CRUST_MANTLE =       475797
 integer, parameter :: NGLOB_OUTER_CORE =        48789
 integer, parameter :: NGLOB_INNER_CORE =         3549

 integer, parameter :: NSPECMAX_ANISO_IC =            1

 integer, parameter :: NSPECMAX_ISO_MANTLE =         7092
 integer, parameter :: NSPECMAX_TISO_MANTLE =         7092
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            1

 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUAT =         7092
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =           45

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =         7092
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =           45

 integer, parameter :: NSPEC_CRUST_MANTLE_STR_AND_ATT =         7092
 integer, parameter :: NSPEC_INNER_CORE_STR_AND_ATT =           45

 integer, parameter :: NSPEC_CRUST_MANTLE_STRAIN_ONLY =         7092
 integer, parameter :: NSPEC_INNER_CORE_STRAIN_ONLY =           45

 integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT =         7092
 integer, parameter :: NSPEC_OUTER_CORE_ADJOINT =          684
 integer, parameter :: NSPEC_INNER_CORE_ADJOINT =           45
 integer, parameter :: NGLOB_CRUST_MANTLE_ADJOINT =       475797
 integer, parameter :: NGLOB_OUTER_CORE_ADJOINT =        48789
 integer, parameter :: NGLOB_INNER_CORE_ADJOINT =         3549
 integer, parameter :: NSPEC_OUTER_CORE_ROT_ADJOINT =          684

 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =         7092
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =          684

 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =       475797

 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.

 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.

 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.

 logical, parameter :: ATTENUATION_VAL = .true.

 logical, parameter :: ATTENUATION_3D_VAL = .false.

 logical, parameter :: ELLIPTICITY_VAL = .true.

 logical, parameter :: GRAVITY_VAL = .true.

 logical, parameter :: OCEANS_VAL = .true.

 logical, parameter :: ROTATION_VAL = .true.
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =          684

 integer, parameter :: NGLOB1D_RADIAL_CM =          189
 integer, parameter :: NGLOB1D_RADIAL_OC =          133
 integer, parameter :: NGLOB1D_RADIAL_IC =           21
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_CM =        10274
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_OC =         2862
 integer, parameter :: NGLOB2DMAX_XMIN_XMAX_IC =          314
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_CM =        10274
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_OC =         2862
 integer, parameter :: NGLOB2DMAX_YMIN_YMAX_IC =          314
 integer, parameter :: NPROC_XI_VAL =           10
 integer, parameter :: NPROC_ETA_VAL =           10
 integer, parameter :: NCHUNKS_VAL =            1
 integer, parameter :: NPROCTOT_VAL =          100
 integer, parameter :: NGLOB2DMAX_XY_VAL =        10274
 integer, parameter :: NUMMSGS_FACES_VAL =           10
 integer, parameter :: NCORNERSCHUNKS_VAL =            1
 integer, parameter :: ATT1 =            5
 integer, parameter :: ATT2 =            5
 integer, parameter :: ATT3 =            5
 integer, parameter :: ATT4 =         7092
 integer, parameter :: ATT5 =           45
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          522
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          522
 integer, parameter :: NSPEC2D_BOTTOM_CM =           36
 integer, parameter :: NSPEC2D_TOP_CM =          576
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =           15
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =           15
 integer, parameter :: NSPEC2D_BOTTOM_IC =            9
 integer, parameter :: NSPEC2D_TOP_IC =            9
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          141
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          144
 integer, parameter :: NSPEC2D_BOTTOM_OC =            9
 integer, parameter :: NSPEC2D_TOP_OC =           36
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1
 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 logical, parameter :: USE_ATTENUATION_MIMIC = .true.
 logical, parameter :: COMPUTE_AND_STORE_STRAIN = .true.
