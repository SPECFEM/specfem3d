!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! read arrays created by the mesher

  subroutine read_arrays_solver(myrank,NSPEC_AB,NGLOB_AB,xstore,ystore,zstore, &
               xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian, &
               flag_sediments,not_fully_in_bedrock,rho_vp,rho_vs,ANISOTROPY, &
               c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
               c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
               c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
               kappastore,mustore,ibool,idoubling,rmass,rmass_ocean_load,LOCAL_PATH,OCEANS)

  implicit none

  include "constants.h"

!  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  integer NSPEC_AB
  integer NGLOB_AB

  logical OCEANS

  character(len=256) LOCAL_PATH

! coordinates in single precision
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian

  logical ANISOTROPY
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO) :: &
            c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
            c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
            c36store,c44store,c45store,c46store,c55store,c56store,c66store

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
  character(len=256) prname

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

! xix
  open(unit=IIN,file=prname(1:len_trim(prname))//'xix.bin',status='old',action='read',form='unformatted')
  read(IIN) xix
  close(IIN)

! xiy
  open(unit=IIN,file=prname(1:len_trim(prname))//'xiy.bin',status='old',action='read',form='unformatted')
  read(IIN) xiy
  close(IIN)

! xiz
  open(unit=IIN,file=prname(1:len_trim(prname))//'xiz.bin',status='old',action='read',form='unformatted')
  read(IIN) xiz
  close(IIN)

! etax
  open(unit=IIN,file=prname(1:len_trim(prname))//'etax.bin',status='old',action='read',form='unformatted')
  read(IIN) etax
  close(IIN)

! etay
  open(unit=IIN,file=prname(1:len_trim(prname))//'etay.bin',status='old',action='read',form='unformatted')
  read(IIN) etay
  close(IIN)

! etaz
  open(unit=IIN,file=prname(1:len_trim(prname))//'etaz.bin',status='old',action='read',form='unformatted')
  read(IIN) etaz
  close(IIN)

! gammax
  open(unit=IIN,file=prname(1:len_trim(prname))//'gammax.bin',status='old',action='read',form='unformatted')
  read(IIN) gammax
  close(IIN)

! gammay
  open(unit=IIN,file=prname(1:len_trim(prname))//'gammay.bin',status='old',action='read',form='unformatted')
  read(IIN) gammay
  close(IIN)

! gammaz
  open(unit=IIN,file=prname(1:len_trim(prname))//'gammaz.bin',status='old',action='read',form='unformatted')
  read(IIN) gammaz
  close(IIN)

! jacobian
  open(unit=IIN,file=prname(1:len_trim(prname))//'jacobian.bin',status='old',action='read',form='unformatted')
  read(IIN) jacobian
  close(IIN)

! read coordinates of the mesh
  open(unit=IIN,file=prname(1:len_trim(prname))//'x.bin',status='old',action='read',form='unformatted')
  read(IIN) xstore
  close(IIN)

  open(unit=IIN,file=prname(1:len_trim(prname))//'y.bin',status='old',action='read',form='unformatted')
  read(IIN) ystore
  close(IIN)

  open(unit=IIN,file=prname(1:len_trim(prname))//'z.bin',status='old',action='read',form='unformatted')
  read(IIN) zstore
  close(IIN)

! ibool
  open(unit=IIN,file=prname(1:len_trim(prname))//'ibool.bin',status='old',action='read',form='unformatted')
  read(IIN) ibool
  close(IIN)

! idoubling
  open(unit=IIN,file=prname(1:len_trim(prname))//'idoubling.bin',status='old',action='read',form='unformatted')
  read(IIN) idoubling
  close(IIN)

! mass matrix
  open(unit=IIN,file=prname(1:len_trim(prname))//'rmass.bin',status='old',action='read',form='unformatted')
  read(IIN) rmass
  close(IIN)

! read additional ocean load mass matrix
  if(OCEANS) then
    open(unit=IIN,file=prname(1:len_trim(prname))//'rmass_ocean_load.bin',status='old',action='read',form='unformatted')
    read(IIN) rmass_ocean_load
    close(IIN)
  endif

! flag_sediments
  open(unit=IIN,file=prname(1:len_trim(prname))//'flag_sediments.bin',status='old',action='read',form='unformatted')
  read(IIN) flag_sediments
  close(IIN)

! not_fully_in_bedrock
  open(unit=IIN,file=prname(1:len_trim(prname))//'not_fully_in_bedrock.bin',status='old',action='read',form='unformatted')
  read(IIN) not_fully_in_bedrock
  close(IIN)

! rho_vs
! Stacey

! rho_vp
  open(unit=IIN,file=prname(1:len_trim(prname))//'rho_vp.bin',status='old',action='read',form='unformatted')
  read(IIN) rho_vp
  close(IIN)

! rho_vs
  open(unit=IIN,file=prname(1:len_trim(prname))//'rho_vs.bin',status='old',action='read',form='unformatted')
  read(IIN) rho_vs
  close(IIN)


! model arrays

! kappa
  open(unit=IIN,file=prname(1:len_trim(prname))//'kappa.bin',status='old',action='read',form='unformatted')
  read(IIN) kappastore
  close(IIN)

! mu
  open(unit=IIN,file=prname(1:len_trim(prname))//'mu.bin',status='old',action='read',form='unformatted')
  read(IIN) mustore
  close(IIN)

  if(ANISOTROPY) then

! c11
     open(unit=IIN,file=prname(1:len_trim(prname))//'c11.bin',status='old',action='read',form='unformatted')
     read(IIN) c11store
     close(IIN)

! c12
     open(unit=IIN,file=prname(1:len_trim(prname))//'c12.bin',status='old',action='read',form='unformatted')
     read(IIN) c12store
     close(IIN)

! c13
     open(unit=IIN,file=prname(1:len_trim(prname))//'c13.bin',status='old',action='read',form='unformatted')
     read(IIN) c13store
     close(IIN)

! c14
     open(unit=IIN,file=prname(1:len_trim(prname))//'c14.bin',status='old',action='read',form='unformatted')
     read(IIN) c14store
     close(IIN)

! c15
     open(unit=IIN,file=prname(1:len_trim(prname))//'c15.bin',status='old',action='read',form='unformatted')
     read(IIN) c15store
     close(IIN)

! c16
     open(unit=IIN,file=prname(1:len_trim(prname))//'c16.bin',status='old',action='read',form='unformatted')
     read(IIN) c16store
     close(IIN)

! c22
     open(unit=IIN,file=prname(1:len_trim(prname))//'c22.bin',status='old',action='read',form='unformatted')
     read(IIN) c22store
     close(IIN)

! c23
     open(unit=IIN,file=prname(1:len_trim(prname))//'c23.bin',status='old',action='read',form='unformatted')
     read(IIN) c23store
     close(IIN)

! c24
     open(unit=IIN,file=prname(1:len_trim(prname))//'c24.bin',status='old',action='read',form='unformatted')
     read(IIN) c24store
     close(IIN)

! c25
     open(unit=IIN,file=prname(1:len_trim(prname))//'c25.bin',status='old',action='read',form='unformatted')
     read(IIN) c25store
     close(IIN)

! c26
     open(unit=IIN,file=prname(1:len_trim(prname))//'c26.bin',status='old',action='read',form='unformatted')
     read(IIN) c26store
     close(IIN)

! c33
     open(unit=IIN,file=prname(1:len_trim(prname))//'c33.bin',status='old',action='read',form='unformatted')
     read(IIN) c33store
     close(IIN)

! c34
     open(unit=IIN,file=prname(1:len_trim(prname))//'c34.bin',status='old',action='read',form='unformatted')
     read(IIN) c34store
     close(IIN)

! c35
     open(unit=IIN,file=prname(1:len_trim(prname))//'c35.bin',status='old',action='read',form='unformatted')
     read(IIN) c35store
     close(IIN)

! c36
     open(unit=IIN,file=prname(1:len_trim(prname))//'c36.bin',status='old',action='read',form='unformatted')
     read(IIN) c36store
     close(IIN)

! c44
     open(unit=IIN,file=prname(1:len_trim(prname))//'c44.bin',status='old',action='read',form='unformatted')
     read(IIN) c44store
     close(IIN)

! c45
     open(unit=IIN,file=prname(1:len_trim(prname))//'c45.bin',status='old',action='read',form='unformatted')
     read(IIN) c45store
     close(IIN)

! c46
     open(unit=IIN,file=prname(1:len_trim(prname))//'c46.bin',status='old',action='read',form='unformatted')
     read(IIN) c46store
     close(IIN)

! c55
     open(unit=IIN,file=prname(1:len_trim(prname))//'c55.bin',status='old',action='read',form='unformatted')
     read(IIN) c55store
     close(IIN)

! c56
     open(unit=IIN,file=prname(1:len_trim(prname))//'c56.bin',status='old',action='read',form='unformatted')
     read(IIN) c56store
     close(IIN)

! c66
     open(unit=IIN,file=prname(1:len_trim(prname))//'c66.bin',status='old',action='read',form='unformatted')
     read(IIN) c66store
     close(IIN)

  endif

  end subroutine read_arrays_solver

