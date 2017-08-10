! Proc  32: Header for mesh information to run static solver
! created by the mesher on 12/11/2015, at 19h 01min

!:::::::::::::::::::: Input parameters :::::::::::::::::::::::::::
!   Background model     :            prem_iso
!   Dominant period [s]  :   10.0000
!   Elements/wavelength  :    1.5000
!   Courant number       :    0.6000
!   Coarsening levels    :         3
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 integer, parameter ::         npol =         4  !            polynomial order
 integer, parameter ::        nelem =      3240  !                   proc. els
 integer, parameter ::       npoint =     81000  !               proc. all pts
 integer, parameter ::    nel_solid =      2712  !             proc. solid els
 integer, parameter ::    nel_fluid =       528  !             proc. fluid els
 integer, parameter :: npoint_solid =     67800  !             proc. solid pts
 integer, parameter :: npoint_fluid =     13200  !             proc. fluid pts
 integer, parameter ::  nglob_fluid =      8702  !            proc. flocal pts
 integer, parameter ::     nel_bdry =        48  ! proc. solid-fluid bndry els
 integer, parameter ::        ndisc =        11  !   # disconts in bkgrd model
 integer, parameter ::   nproc_mesh =        32  !        number of processors
 integer, parameter :: lfbkgrdmodel =         8  !   length of bkgrdmodel name

!:::::::::::::::::::: Output parameters ::::::::::::::::::::::::::
!   Time step [s]        :    0.0733
!   Min(h/vp),dt/courant :    0.6693    0.4890
!   max(h/vs),T0/wvlngth :    6.6391    6.6667
!   Inner core r_min [km]: 1159.7965
!   Max(h) r/ns(icb) [km]:   14.9901
!   Max(h) precalc.  [km]:   17.8430
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

