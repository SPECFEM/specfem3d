! Proc  32: Header for mesh information to run static solver
! created by the mesher on 06/11/2015, at 18h 27min
 
!:::::::::::::::::::: Input parameters :::::::::::::::::::::::::::
!   Background model     :              iasp91
!   Inner-core shear wave:         T
!   Dominant period [s]  :    5.0000
!   Elements/wavelength  :    1.5000
!   Courant number       :    0.6000
!   Coarsening levels    :         3
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
 integer, parameter ::         npol =         4  !            polynomial order
 integer, parameter ::        nelem =     12124  !                   proc. els
 integer, parameter ::       npoint =    303100  !               proc. all pts
 integer, parameter ::    nel_solid =     10164  !             proc. solid els
 integer, parameter ::    nel_fluid =      1960  !             proc. fluid els
 integer, parameter :: npoint_solid =    254100  !             proc. solid pts
 integer, parameter :: npoint_fluid =     49000  !             proc. fluid pts
 integer, parameter ::  nglob_fluid =     31794  !            proc. flocal pts
 integer, parameter ::     nel_bdry =        84  ! proc. solid-fluid bndry els
 integer, parameter ::        ndisc =        11  !   # disconts in bkgrd model
 integer, parameter ::   nproc_mesh =        32  !        number of processors
 integer, parameter :: lfbkgrdmodel =         6  !   length of bkgrdmodel name
 
!:::::::::::::::::::: Output parameters ::::::::::::::::::::::::::
!   Time step [s]        :    0.0381
!   Min(h/vp),dt/courant :    0.3657    0.2538
!   max(h/vs),T0/wvlngth :    3.3274    3.3333
!   Inner core r_min [km]:  989.2201
!   Max(h) r/ns(icb) [km]:    8.5342
!   Max(h) precalc.  [km]:    8.7542
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
