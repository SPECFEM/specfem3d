!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 1
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology October 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! read arrays created by the mesher

  subroutine read_arrays_solver(myrank,xstore,ystore,zstore, &
               xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian, &
               flag_sediments,not_fully_in_bedrock,rho_vp,rho_vs, &
               kappastore,mustore,ibool,idoubling,rmass,rmass_ocean_load,LOCAL_PATH,OCEANS)

  implicit none

  include "constants.h"

  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  logical OCEANS

  character(len=150) LOCAL_PATH

! coordinates in single precision
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian

! material properties
  real(kind=CUSTOM_REAL) kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL) mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB)

! flag for sediments
  logical not_fully_in_bedrock(NSPEC_AB)
  logical flag_sediments(NGLLX,NGLLY,NGLLZ,NSPEC_AB)

! Stacey
  real(kind=CUSTOM_REAL) rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL) rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB)

! mass matrix and additional ocean load mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: rmass,rmass_ocean_load

! global addressing
  integer ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB)

  integer idoubling(NSPEC_AB)

! processor identification
  character(len=150) prname

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

! xix
  open(unit=IIN,file=prname(1:len_trim(prname))//'xix.bin',status='old',form='unformatted')
  read(IIN) xix
  close(IIN)

! xiy
  open(unit=IIN,file=prname(1:len_trim(prname))//'xiy.bin',status='old',form='unformatted')
  read(IIN) xiy
  close(IIN)

! xiz
  open(unit=IIN,file=prname(1:len_trim(prname))//'xiz.bin',status='old',form='unformatted')
  read(IIN) xiz
  close(IIN)

! etax
  open(unit=IIN,file=prname(1:len_trim(prname))//'etax.bin',status='old',form='unformatted')
  read(IIN) etax
  close(IIN)

! etay
  open(unit=IIN,file=prname(1:len_trim(prname))//'etay.bin',status='old',form='unformatted')
  read(IIN) etay
  close(IIN)

! etaz
  open(unit=IIN,file=prname(1:len_trim(prname))//'etaz.bin',status='old',form='unformatted')
  read(IIN) etaz
  close(IIN)

! gammax
  open(unit=IIN,file=prname(1:len_trim(prname))//'gammax.bin',status='old',form='unformatted')
  read(IIN) gammax
  close(IIN)

! gammay
  open(unit=IIN,file=prname(1:len_trim(prname))//'gammay.bin',status='old',form='unformatted')
  read(IIN) gammay
  close(IIN)

! gammaz
  open(unit=IIN,file=prname(1:len_trim(prname))//'gammaz.bin',status='old',form='unformatted')
  read(IIN) gammaz
  close(IIN)

! jacobian
  open(unit=IIN,file=prname(1:len_trim(prname))//'jacobian.bin',status='old',form='unformatted')
  read(IIN) jacobian
  close(IIN)

! read coordinates of the mesh
  open(unit=IIN,file=prname(1:len_trim(prname))//'x.bin',status='old',form='unformatted')
  read(IIN) xstore
  close(IIN)

  open(unit=IIN,file=prname(1:len_trim(prname))//'y.bin',status='old',form='unformatted')
  read(IIN) ystore
  close(IIN)

  open(unit=IIN,file=prname(1:len_trim(prname))//'z.bin',status='old',form='unformatted')
  read(IIN) zstore
  close(IIN)

! ibool
  open(unit=IIN,file=prname(1:len_trim(prname))//'ibool.bin',status='old',form='unformatted')
  read(IIN) ibool
  close(IIN)

! idoubling
  open(unit=IIN,file=prname(1:len_trim(prname))//'idoubling.bin',status='old',form='unformatted')
  read(IIN) idoubling
  close(IIN)

! mass matrix
  open(unit=IIN,file=prname(1:len_trim(prname))//'rmass.bin',status='old',form='unformatted')
  read(IIN) rmass
  close(IIN)

! read additional ocean load mass matrix
  if(OCEANS) then
    open(unit=IIN,file=prname(1:len_trim(prname))//'rmass_ocean_load.bin',status='old',form='unformatted')
    read(IIN) rmass_ocean_load
    close(IIN)
  endif

! flag_sediments
  open(unit=IIN,file=prname(1:len_trim(prname))//'flag_sediments.bin',status='old',form='unformatted')
  read(IIN) flag_sediments
  close(IIN)

! not_fully_in_bedrock
  open(unit=IIN,file=prname(1:len_trim(prname))//'not_fully_in_bedrock.bin',status='old',form='unformatted')
  read(IIN) not_fully_in_bedrock
  close(IIN)

! rho_vs
! Stacey

! rho_vp
  open(unit=IIN,file=prname(1:len_trim(prname))//'rho_vp.bin',status='old',form='unformatted')
  read(IIN) rho_vp
  close(IIN)

! rho_vs
  open(unit=IIN,file=prname(1:len_trim(prname))//'rho_vs.bin',status='old',form='unformatted')
  read(IIN) rho_vs
  close(IIN)


! model arrays

! kappa
  open(unit=IIN,file=prname(1:len_trim(prname))//'kappa.bin',status='old',form='unformatted')
  read(IIN) kappastore
  close(IIN)

! mu
  open(unit=IIN,file=prname(1:len_trim(prname))//'mu.bin',status='old',form='unformatted')
  read(IIN) mustore
  close(IIN)

  end subroutine read_arrays_solver

