! Proc  32: Header for mesh information to run static solver
! created by the mesher on 09/14/2015, at 18h 17min

!:::::::::::::::::::: Input parameters :::::::::::::::::::::::::::
!   Background model     :            prem_iso
!   Dominant period [s]  :    5.0000
!   Elements/wavelength  :    1.5000
!   Courant number       :    0.6000
!   Coarsening levels    :         3
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 integer, parameter ::         npol =         4  !            polynomial order
 integer, parameter ::        nelem =     11858  !                   proc. els
 integer, parameter ::       npoint =    296450  !               proc. all pts
 integer, parameter ::    nel_solid =      9912  !             proc. solid els
 integer, parameter ::    nel_fluid =      1946  !             proc. fluid els
 integer, parameter :: npoint_solid =    247800  !             proc. solid pts
 integer, parameter :: npoint_fluid =     48650  !             proc. fluid pts
 integer, parameter ::  nglob_fluid =     31574  !            proc. flocal pts
 integer, parameter ::     nel_bdry =        84  ! proc. solid-fluid bndry els
 integer, parameter ::        ndisc =        11  !   # disconts in bkgrd model
 integer, parameter ::   nproc_mesh =        32  !        number of processors
 integer, parameter :: lfbkgrdmodel =         8  !   length of bkgrdmodel name

!:::::::::::::::::::: Output parameters ::::::::::::::::::::::::::
!   Time step [s]        :    0.0389
!   Min(h/vp),dt/courant :    0.3739    0.2596
!   max(h/vs),T0/wvlngth :    3.3272    3.3333
!   Inner core r_min [km]: 1008.1308
!   Max(h) r/ns(icb) [km]:    8.5657
!   Max(h) precalc.  [km]:    8.9215
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

