!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

!--------------------------------------------------------------------------------------------------
!
! PREM [Dziewonski and Anderson, 1981].
!
! A. M. Dziewonski and D. L. Anderson.
! Preliminary reference Earth model.
! Phys. Earth Planet. Inter., 25:297–356, 1981.
!
! Isotropic (iso) and transversely isotropic (aniso) version of the
! spherically symmetric Preliminary Reference Earth Model
!
!--------------------------------------------------------------------------------------------------


  subroutine model_1D_prem_iso(xmesh,ymesh,zmesh,rho_prem,vp_prem,vs_prem,qmu_atten)

!
! isotropic prem model
!

! given a GLL point, returns super-imposed velocity model values

  use generate_databases_par,only: nspec => NSPEC_AB,ibool
  use create_regions_mesh_ext_par
  implicit none

  ! GLL point location
  double precision, intent(in) :: xmesh,ymesh,zmesh

  ! material properties
  real(kind=CUSTOM_REAL), intent(inout) :: rho_prem,vp_prem,vs_prem,qmu_atten

  ! local parameters
  real(kind=CUSTOM_REAL) :: xloc,yloc,zloc
  real(kind=CUSTOM_REAL) :: depth
  real(kind=CUSTOM_REAL) :: elevation,distmin
  double precision :: x,rho,drhodr,vp,vs,Qkappa,Qmu
  double precision :: &
      R_EARTH_M,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R771,R600,R670,R400,R220,R80,RMOHO,RMIDDLE_CRUST,ROCEAN
  double precision :: r

  ! uses crustal values from other models (like crust2.0) than prem
  ! set to .false. to use PREM crustal values, otherwise will take mantle values up to surface
  logical,parameter :: CRUSTAL = .false.
  ! avoids crustal values, uses Moho values up to the surface
  logical,parameter :: SUPPRESS_CRUSTAL_MESH = .false.
  ! same properties everywhere in PREM crust if we decide to define only one layer in the crust
  logical,parameter :: ONE_CRUST = .false.

  ! GLL point location converted to real
  xloc = xmesh
  yloc = ymesh
  zloc = zmesh

  ! get approximate topography elevation at target coordinates
  distmin = HUGEVAL
  elevation = 0.0
  call get_topo_elevation_free_closest(xloc,yloc,elevation,distmin, &
                    nspec,nglob_dummy,ibool,xstore_dummy,ystore_dummy,zstore_dummy, &
                    num_free_surface_faces,free_surface_ispec,free_surface_ijk)

  ! depth in Z-direction
  if( distmin < HUGEVAL ) then
    depth = elevation - zloc
  else
    depth = - zloc
  endif

  ! PREM layers (in m)
  R_EARTH_M = 6371000.d0
  ROCEAN = 6368000.d0
  RMIDDLE_CRUST = 6356000.d0
  RMOHO = 6346600.d0
  R80  = 6291000.d0
  R220 = 6151000.d0
  R400 = 5971000.d0
  R600 = 5771000.d0
  R670 = 5701000.d0
  R771 = 5600000.d0
  RTOPDDOUBLEPRIME = 3630000.d0
  RCMB = 3480000.d0
  RICB = 1221000.d0

  ! compute real physical radius in meters
  r = R_EARTH - dble(depth)

  ! normalized radius
  x = r / R_EARTH

  ! given a normalized radius x, gives the non-dimensionalized density rho,
  ! speeds vp and vs, and the quality factors Qkappa and Qmu

  ! initializes
  rho = 0.d0
  vp = 0.d0
  vs = 0.d0
  Qmu = 0.d0
  Qkappa = 0.d0
  drhodr = 0.d0

  !
  !--- inner core
  !
  if(r >= 0.d0 .and. r <= RICB) then
    drhodr=-2.0d0*8.8381d0*x
    rho=13.0885d0-8.8381d0*x*x
    vp=11.2622d0-6.3640d0*x*x
    vs=3.6678d0-4.4475d0*x*x
    Qmu=84.6d0
    Qkappa=1327.7d0
  !
  !--- outer core
  !
  else if(r > RICB .and. r <= RCMB) then
    drhodr=-1.2638d0-2.0d0*3.6426d0*x-3.0d0*5.5281d0*x*x
    rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
    vp=11.0487d0-4.0362d0*x+4.8023d0*x*x-13.5732d0*x*x*x
    vs=0.0d0
    Qmu=0.0d0
    Qkappa=57827.0d0
  !
  !--- D" at the base of the mantle
  !
  else if(r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
    drhodr=-6.4761d0+2.0d0*5.5283d0*x-3.0d0*3.0807d0*x*x
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=15.3891d0-5.3181d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vs=6.9254d0+1.4672d0*x-2.0834d0*x*x+0.9783d0*x*x*x
    Qmu=312.0d0
    Qkappa=57827.0d0
  !
  !--- mantle: from top of D" to d670
  !
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
    drhodr=-6.4761d0+2.0d0*5.5283d0*x-3.0d0*3.0807d0*x*x
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=24.9520d0-40.4673d0*x+51.4832d0*x*x-26.6419d0*x*x*x
    vs=11.1671d0-13.7818d0*x+17.4575d0*x*x-9.2777d0*x*x*x
    Qmu=312.0d0
    Qkappa=57827.0d0
  else if(r > R771 .and. r <= R670) then
    drhodr=-6.4761d0+2.0d0*5.5283d0*x-3.0d0*3.0807d0*x*x
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=29.2766d0-23.6027d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vs=22.3459d0-17.2473d0*x-2.0834d0*x*x+0.9783d0*x*x*x
    Qmu=312.0d0
    Qkappa=57827.0d0
  !
  !--- mantle: above d670
  !
  else if(r > R670 .and. r <= R600) then
    drhodr=-1.4836d0
    rho=5.3197d0-1.4836d0*x
    vp=19.0957d0-9.8672d0*x
    vs=9.9839d0-4.9324d0*x
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R600 .and. r <= R400) then
    drhodr=-8.0298d0
    rho=11.2494d0-8.0298d0*x
    vp=39.7027d0-32.6166d0*x
    vs=22.3512d0-18.5856d0*x
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R400 .and. r <= R220) then
    drhodr=-3.8045d0
    rho=7.1089d0-3.8045d0*x
    vp=20.3926d0-12.2569d0*x
    vs=8.9496d0-4.4597d0*x
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R220 .and. r <= R80) then
    drhodr=0.6924d0
    rho=2.6910d0+0.6924d0*x
    vp=4.1875d0+3.9382d0*x
    vs=2.1519d0+2.3481d0*x
    Qmu=80.0d0
    Qkappa=57827.0d0
  else
  if(CRUSTAL .and. .not. SUPPRESS_CRUSTAL_MESH) then
! fill with PREM mantle and later add CRUST2.0
    if(r > R80) then
      ! density/velocity from mantle just below moho
      drhodr=0.6924d0
      rho=2.6910d0+0.6924d0*x
      vp=4.1875d0+3.9382d0*x
      vs=2.1519d0+2.3481d0*x
      ! shear attenuation for R80 to surface
      Qmu=600.0d0
      Qkappa=57827.0d0
    endif
  else
! use PREM crust
    if(r > R80 .and. r <= RMOHO) then
      drhodr=0.6924d0
      rho=2.6910d0+0.6924d0*x
      vp=4.1875d0+3.9382d0*x
      vs=2.1519d0+2.3481d0*x
      Qmu=600.0d0
      Qkappa=57827.0d0

    else if (SUPPRESS_CRUSTAL_MESH) then
      !! DK DK extend the Moho up to the surface instead of the crust
      drhodr=0.6924d0
      rho = 2.6910d0+0.6924d0*(RMOHO / R_EARTH)
      vp = 4.1875d0+3.9382d0*(RMOHO / R_EARTH)
      vs = 2.1519d0+2.3481d0*(RMOHO / R_EARTH)
      Qmu=600.0d0
      Qkappa=57827.0d0

    else if(r > RMOHO .and. r <= RMIDDLE_CRUST) then
      drhodr=0.0d0
      rho=2.9d0
      vp=6.8d0
      vs=3.9d0
      Qmu=600.0d0
      Qkappa=57827.0d0

      ! same properties everywhere in PREM crust if we decide to define only one layer in the crust
      if(ONE_CRUST) then
        drhodr=0.0d0
        rho=2.6d0
        vp=5.8d0
        vs=3.2d0
        Qmu=600.0d0
        Qkappa=57827.0d0
      endif

    else if(r > RMIDDLE_CRUST .and. r <= ROCEAN) then
      drhodr=0.0d0
      rho=2.6d0
      vp=5.8d0
      vs=3.2d0
      Qmu=600.0d0
      Qkappa=57827.0d0
    ! for density profile for gravity, we do not check that r <= R_EARTH
    else if(r > ROCEAN) then
      drhodr=0.0d0
      rho=2.6d0
      vp=5.8d0
      vs=3.2d0
      Qmu=600.0d0
      Qkappa=57827.0d0

    endif
  endif
  endif

  ! scales values to SI units ( m/s, kg/m**3)
  rho_prem=rho*1000.0d0
  vp_prem=vp*1000.0d0
  vs_prem=vs*1000.0d0

  qmu_atten = Qmu

  end subroutine model_1D_prem_iso


!
!-------------------------------------------------------------------------------------------------
!

!PB ROUTINE WHICH ASSINGS AT EACH GLL POINT THE VALUES OF v_p v_s AND rho TAKEN FROM THE PREM MODEL
!PB IT'S CALLED IN get_model.f90(:149) WITHIN THE generate_databases PROCESS

  subroutine model_1D_PREM_routine_PB(xloc,yloc,zloc,ro_prem,vp_prem,vs_prem,idom)

  use generate_databases_par,only: CUSTOM_REAL
  implicit none

  double precision, intent(in) :: xloc,yloc,zloc
  real(kind=CUSTOM_REAL),intent(inout) :: ro_prem,vp_prem,vs_prem
  integer, intent(in)          :: idom

  ! local parameters
  double precision             :: r0,r,x_prem
  !double precision             :: ro_prem,vp_prem,vs_prem
  !character(len=3), intent(in) :: param !rho, vs,vp

  r0=sqrt(xloc**2+yloc**2+zloc**2)
  !  print*,'xloc,yloc,zloc,r0,idom',xloc,yloc,zloc,r0,idom
  r=r0/1000.

  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )
  !  print*,'x_prem',x_prem
  IF(idom==1)THEN        ! upper crustal layer
     !print*,'I am in domain 1'
     ro_prem=2.6
     vp_prem=5.8
     vs_prem=3.2
  !    print*,'ro,vp,vs form domain 1',ro_prem,vp_prem,vs_prem
  else if(idom==2)THEN
     ro_prem=2.9                       ! lower crustal layer
     vp_prem=6.8
     vs_prem=3.9
  else if(idom==3)THEN
     ro_prem=2.691+.6924*x_prem             ! upper mantle
     vp_prem=4.1875+3.9382*x_prem
     vs_prem=2.1519+2.3481*x_prem
  else if(idom==4)THEN
     ro_prem=7.1089-3.8045*x_prem
     vp_prem=20.3926-12.2569*x_prem
     vs_prem=8.9496-4.4597*x_prem
  else if(idom==5)THEN
     ro_prem=11.2494-8.0298*x_prem
     vp_prem=39.7027-32.6166*x_prem
     vs_prem=22.3512-18.5856*x_prem
  else if(idom==6)THEN
     ro_prem=5.3197-1.4836*x_prem
     vp_prem=19.0957-9.8672*x_prem
     vs_prem=9.9839-4.9324*x_prem
  else if(idom==7)THEN   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
  else if(idom==8)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vs_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
  else if(idom==9)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+.9783*x_prem**3
  else if(idom==10)THEN  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vp_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vs_prem=0.00
  else if(idom==11)THEN                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vp_prem=11.2622-6.3640*x_prem**2
     vs_prem=3.6678-4.4475*x_prem**2
  ENDIF

  !    print*,'ro,vp,vs passed from routine and not multiplied',ro_prem,vp_prem,vs_prem

  ro_prem=ro_prem*1000
  vp_prem=vp_prem*1000
  vs_prem=vs_prem*1000

  !   print*,'ro,vp,vs passed from PREM_ROUTINE multiplied by 1000',ro_prem,vp_prem,vs_prem

  !  if (param=='rho') then
  !     prem_sub=ro_prem*1000.
  !  else if (param=='v_p') then
  !     prem_sub=vp_prem*1000.
  !  else if (param=='v_s') then
  !     prem_sub=vs_prem*1000.
  !  else
  !     write(6,*)'ERROR IN PREM_SUB FUNCTION:',param,'NOT AN OPTION'
  !     stop

  end subroutine model_1D_PREM_routine_PB

