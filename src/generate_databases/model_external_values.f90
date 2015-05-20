!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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
! generic model file
!
! note: the idea is to super-impose velocity model values on the GLL points,
!          additional to the ones assigned on the CUBIT mesh
!
! most of the routines here are place-holders, please add/implement your own routines
!
!--------------------------------------------------------------------------------------------------


  module external_model

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! only here to illustrate an example
  !  type model_external_variables
  !    sequence
  !    double precision dvs(0:dummy_size)
  !  end type model_external_variables
  !  type (model_external_variables) MEXT_V

  ! VM VM my model for DSM coupling
  use constants
   ! VM VM
  double precision, dimension (:,:), allocatable :: vpv_1D,vsv_1D,density_1D
  double precision, dimension (:), allocatable :: zlayer
  double precision, dimension (:), allocatable :: smooth_vp,smooth_vs
  integer :: ilayer,nlayer,ncoeff,ndeg_poly
  double precision :: ZREF,OLON,OLAT

  end module external_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_broadcast(myrank)

! standard routine to setup model

  use external_model

  use constants

  use generate_databases_par, only: COUPLE_WITH_EXTERNAL_CODE

  implicit none

  integer :: myrank

  ! local parameters
  integer :: idummy

  ! dummy to ignore compiler warnings
  idummy = myrank

!---
!
! ADD YOUR MODEL HERE
!
  if (COUPLE_WITH_EXTERNAL_CODE) call read_external_model_coupling()

!---

  ! the variables read are declared and stored in structure MEXT_V
  !if (myrank == 0) call read_external_model()

  ! broadcast the information read on the master to the nodes
  !call bcast_all_dp(MEXT_V%dvs, size(MEXT_V%dvs))

  end subroutine model_external_broadcast

!
!-------------------------------------------------------------------------------------------------
!
!
!  subroutine read_external_model()
!
!  use external_model
!
!  use constants
!
!  implicit none
!
!---
!
! ADD YOUR MODEL HERE
!
!---
!
!  end subroutine read_external_model
!
!-------------------------------------------------------------------------------------------------
!
!
  subroutine read_external_model_coupling()
!
  use external_model !! VM VM custom subroutine for coupling with DSM
!
   implicit none
!
!---
!
! ADD YOUR MODEL HERE
!
    character(len=256):: filename
    integer i,cc
    double precision aa,bb

    filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'/coeff_poly_deg12'
    open(27,file=trim(filename))
    read(27,*) ndeg_poly
    allocate(smooth_vp(0:ndeg_poly),smooth_vs(0:ndeg_poly))
    do i=ndeg_poly,0,-1
       read(27,*) aa,bb,cc
       smooth_vp(i) = aa
       smooth_vs(i) = bb
       !write(*,*) a,b
    enddo

    filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'/model_1D.in'
    open(27,file=trim(filename))
    read(27,*) nlayer,ncoeff
    allocate(vpv_1D(nlayer,ncoeff))
    allocate(vsv_1D(nlayer,ncoeff))
    allocate(density_1D(nlayer,ncoeff))
    allocate(zlayer(nlayer))
    do i=1,nlayer
       read(27,*) zlayer(i)
       read(27,*) vpv_1D(i,:)
       read(27,*) vsv_1D(i,:)
       read(27,*) density_1D(i,:)
    enddo
    read(27,*) ZREF
    read(27,*) OLON,OLAT

!---
!
  end subroutine read_external_model_coupling
!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_values(xmesh,ymesh,zmesh,rho,vp,vs,qkappa_atten,qmu_atten,iflag_aniso,idomain_id )

! given a GLL point, returns super-imposed velocity model values

  use generate_databases_par,only: nspec => NSPEC_AB,ibool,HUGEVAL,TINYVAL,IDOMAIN_ELASTIC,COUPLE_WITH_EXTERNAL_CODE

  use create_regions_mesh_ext_par

!  use external_model

  implicit none

  ! GLL point
  double precision, intent(in) :: xmesh,ymesh,zmesh

  ! density, Vp and Vs
  real(kind=CUSTOM_REAL) :: vp,vs,rho

  ! attenuation flag
  real(kind=CUSTOM_REAL) :: qkappa_atten,qmu_atten

  ! anisotropy flag
  integer :: iflag_aniso

  ! acoustic/elastic/.. domain flag ( 1 = acoustic / 2 = elastic / ... )
  integer :: idomain_id

  ! local parameters
  double precision :: radius, x,y,z
  real(kind=CUSTOM_REAL) :: xmin,xmax,ymin,ymax,zmin,zmax,x_target,y_target
  real(kind=CUSTOM_REAL) :: depth
  real(kind=CUSTOM_REAL) :: elevation,distmin

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! GLL point location converted to real
  x = xmesh
  y = ymesh
  z = zmesh

  if (COUPLE_WITH_EXTERNAL_CODE) then
    call  model_1D(x,y,z,rho,vp,vs,radius)
  else
    ! note: z coordinate will be negative below surface
    !          convention is z-axis points up

    ! model dimensions
    xmin = 0._CUSTOM_REAL ! minval(xstore_dummy)
    xmax = 134000._CUSTOM_REAL ! maxval(xstore_dummy)
    ymin = 0._CUSTOM_REAL  !minval(ystore_dummy)
    ymax = 134000._CUSTOM_REAL ! maxval(ystore_dummy)
    zmin = 0._CUSTOM_REAL ! minval(zstore_dummy)
    zmax = 60000._CUSTOM_REAL ! maxval(zstore_dummy)
    x_target = x
    y_target = y

    ! get approximate topography elevation at target coordinates from free surface
    call get_topo_elevation_free_closest(x_target,y_target,elevation,distmin, &
         nspec,nglob_dummy,ibool,xstore_dummy,ystore_dummy,zstore_dummy, &
         num_free_surface_faces,free_surface_ispec,free_surface_ijk)


    ! depth in Z-direction
    if (distmin < HUGEVAL) then
       depth = elevation - z
    else
       depth = zmin - z
    endif

    ! normalizes depth between 0 and 1
    if (abs( zmax - zmin ) > TINYVAL) depth = depth / (zmax - zmin)

    ! initial values (in m/s and kg/m^3)
    rho = 2691.0_CUSTOM_REAL
    vp = 4187.5_CUSTOM_REAL
    vs = 2151.9_CUSTOM_REAL

    ! adds a velocity depth gradient
    ! (e.g. from PREM mantle gradients:
    !     vp : 3.9382*6371/5.5
    !     vs : 2.3481*6371/5.5
    !     rho : 0.6924*6371/5.5 )
    rho = rho + 802._CUSTOM_REAL * depth
    vp = vp + 4562._CUSTOM_REAL * depth
    vs = vs + 2720._CUSTOM_REAL * depth

    ! attenuation: PREM crust value
    qmu_atten=600._CUSTOM_REAL

    ! no Q_kappa in this model
    qkappa_atten = 9999._CUSTOM_REAL

    ! no anisotropy
    iflag_aniso = 0

    ! elastic material
    idomain_id = IDOMAIN_ELASTIC

  endif

  end subroutine model_external_values

!----------------------------------------------------------------
!! !! ================= VM VM CUSTOM SUBROUTINE FOR DSM COUPLING
!----------------------------------------------------------------

  subroutine FindLayer(x,y,z)

    use external_model
    implicit none
    integer il
    double precision radius
    double precision :: x,y,z
    radius =  dsqrt(x**2 + y**2 + (z+zref)**2) / 1000.d0

    !write(124,*) 'RADIUS ',radius,x,y,z,z+zref,zref
    il = 1
    do while (radius > zlayer(il).and.il<nlayer)
       il = il + 1
    enddo
    il = il - 1
    ilayer = il

    !write(124,*) 'r, i : ',z,zref,radius, ilayer
  end subroutine FindLayer

!----------------------------------------------------------------

 subroutine model_1D(x_eval,y_eval,z_eval, &
                             rho_final,vp_final,vs_final,r1)
    use external_model
    implicit none
    double precision r1,radius
    double precision rho,vp,vs
    double precision x_eval,y_eval,z_eval
    real(kind=CUSTOM_REAL) rho_final,vp_final,vs_final
    double precision Interpol,Xtol


    Xtol=1d-2

    radius = dsqrt(x_eval**2 + y_eval**2 + (z_eval+zref)**2)
    radius = radius / 1000.d0
    r1=radius

    ! get vp,vs and rho
    radius = radius / zlayer(nlayer)
    vp = Interpol(vpv_1D,ilayer,radius,nlayer)
    vs = Interpol(vsv_1D,ilayer,radius,nlayer)
    rho = Interpol(density_1D,ilayer,radius,nlayer)

    vp_final = vp * 1000.d0
    vs_final = vs * 1000.d0
    rho_final = rho * 1000.d0

  end subroutine model_1D

!----------------------------------------------------------------

  function Interpol(v,i,x,nl)
    implicit none
    integer i,nl
    double precision Interpol,x,v(nl,4)
    Interpol = v(i,1)+x*(v(i,2)+x*(v(i,3)+x*v(i,4)))
  end function Interpol

!----------------------------------------------------------------
!! CD CD ================= BELOW : VM VM SUBROUTINES FOR ??
!----------------------------------------------------------------

  subroutine model_1D_smooth(x_eval,y_eval,z_eval, &
                             rho_final,vp_final,vs_final,r1)

    use external_model
    implicit none
    double precision r1,radius,radius_earth_km
    double precision rho,vp,vs
    double precision x_eval,y_eval,z_eval
    real(kind=CUSTOM_REAL) rho_final,vp_final,vs_final
    double precision Interpol,Xtol
    double precision vp1,vs1,rho1
    double precision depth
    integer ideg

!!$   !! HARDCODED ---------------
!!$    !! ZONE 1 externe
!!$    xmin_1=-108000.d0
!!$    xmax_1= 108000.d0
!!$
!!$    ymin_1=-108000.d0
!!$    ymax_1= 108000.d0
!!$
!!$    zmin_1 = 5000.d0
!!$
!!$    !! ZONE 2 intermediaire
!!$    xmin_2= -85000.d0
!!$    xmax_2=  85000.d0
!!$
!!$    ymin_2= -85000.d0
!!$    ymax_2=  85000.d0
!!$
!!$    zmin_2 = 10000.d0



    radius_earth_km = 6371.d0
    Xtol=1d-2

    radius = dsqrt(x_eval**2 + y_eval**2 + (z_eval+zref)**2)
    radius = radius / 1000.d0
    r1=radius

    depth =  radius_earth_km - r1

    vp1=0.d0;vs1=0.d0
    do ideg=0,ndeg_poly
      vp1 = vp1 + smooth_vp(ideg) * (depth  * 1000.d0 )** ideg
      vs1 = vs1 + smooth_vs(ideg) * (depth  * 1000.d0 )** ideg
    end do

    vp1=1.d0/vp1
    vs1=1.d0/vs1
    rho1=vp1/2.34 ! rho_final/vp_final

    ! get vp,vs and rho
    !radius = radius / zlayer(nlayer)
    !vp = Interpol(vpv_1D,ilayer,radius,nlayer)
    !vs = Interpol(vsv_1D,ilayer,radius,nlayer)
    !rho = Interpol(density_1D,ilayer,radius,nlayer)

    !vp_final = vp * 1000.d0
    !vs_final = vs * 1000.d0
    !rho_final = rho * 1000.d0
!
   call tapper_3D( x_eval, y_eval, z_eval, vp_final, vs_final, rho_final, vp1, vs1, rho1)


  end subroutine model_1D_smooth



  subroutine model_3d_vm(x_eval,y_eval,z_eval, &
                             rho_final,vp_final,vs_final)
    
    use external_model
    implicit none
    double precision r1,radius
    double precision rho,vp,vs
    double precision x_eval,y_eval,z_eval
    double precision x,y,z
    real(kind=CUSTOM_REAL) rho_final,vp_final,vs_final
    double precision Interpol,Xtol
    integer idepth,ilon,ilat,izone,iflag
    double precision  :: tapperx, tappery, tapperz, tapper
    double precision :: vpi,vsi,rhoi, ri, rj, rk
!!$    !! HARDCODED ---------------
!!$    !! ZONE 1 externe
!!$    xmin_1=-108000.d0
!!$    xmax_1= 108000.d0
!!$
!!$    ymin_1=-108000.d0
!!$    ymax_1= 108000.d0
!!$    
!!$    zmin_1 = 5000.d0
!!$    
!!$    !! ZONE 2 intermediaire
!!$    xmin_2= -85000.d0
!!$    xmax_2=  85000.d0
!!$
!!$    ymin_2= -85000.d0
!!$    ymax_2=  85000.d0
!!$    
!!$    zmin_2 = 10000.d0

    x=x_eval
    y=y_eval
    z=z_eval
    Xtol=1d-2


    radius = dsqrt(x_eval**2 + y_eval**2 + (z_eval+zref)**2)
    idepth = -floor((radius - (Radius_Earth - z3d(1,1,1) )) / hz3d) + 1  
    ilon   = floor( (x_eval - x3d(1,1,1)) / hx3d ) + 1
    ilat  =  floor( (y_eval - y3d(1,1,1)) / hy3d ) + 1
    

   
    
    !write(*,*) 

    if (idepth < 1) idepth = 1
    if (idepth > nz3d) idepth = nz3d

    if (ilon < 1) ilon = 1
    if (ilon > nx3d) ilon = nx3d

    if (ilat < 1) ilat = 1
    if (ilat > nx3d) ilat = ny3d


    if (idepth > 1 .and. idepth < nz3d .and. ilon < nx3d .and. ilat < ny3d) then
       ri = (x_eval - x3d(ilon, ilat, idepth)) /  hx3d
       rj = (y_eval - y3d(ilon, ilat, idepth)) /  hy3d
       rk = -(radius -  (Radius_Earth - z3d(ilon,ilat,idepth-1) )) / hz3d  
    !write(*,*) ri,rj,rk
    !! tri-linear interpolation
       call tri_interp(ilon,ilat,idepth-1,ri,rj,rk,vpi,vsi,rhoi)
       !if (vpi==0.) write(*,*) ri,rj,rk
    else
       vpi=vp3d(ilon,ilat,idepth)
       vsi=vs3d(ilon,ilat,idepth)
       rhoi=rho3d(ilon,ilat,idepth)
    end if

    call tapper_3D( x, y, z, vp_final, vs_final, rho_final, vpi, vsi, rhoi)

!!!!!  tapper in side --
!    call in_which_zone_I_am(x,y,z,izone) 
!
!    if (izone == 1) then
!       vp_final  = vpi  ! vp3d(ilon,ilat,idepth)   ! * 1000.d0
       !write(*,*) idepth, vp_final
!       vs_final  = vsi  ! vs3d(ilon,ilat,idepth)   ! * 1000.d0
!       rho_final = rhoi !rho3d(ilon,ilat,idepth) ! * 1000.d0
!    end if
!     if (izone == 2) then
!       tapperx=1.d0
!       tappery=1.d0
!       tapperz=1.d0
!       tapper=1.d0
       ! add tapper
!       if     (x > xmin_1 .and. x < xmin_2) then
!          iflag = 1
!          call tapper_1D(x,xmin_1,xmin_2,tapperx, iflag)
          

!       elseif (x < xmax_1 .and. x > xmax_2) then
!          iflag = -1
!          call tapper_1D(x,xmax_1,xmax_2,tapperx, iflag)

!       end if
       
!       if     (y > ymin_1 .and. y < ymin_2) then
!          iflag = 1
!          call tapper_1D(y,ymin_1,ymin_2,tappery, iflag)
          
!       elseif (y < ymax_1 .and. y > ymax_2) then
!          iflag = -1
!          call tapper_1D(y,ymax_1,ymax_2,tappery, iflag)
!       end if

!       if     (z > zmin_1 .and. z < zmin_2) then
!          iflag = 1
!          call tapper_1D(z,zmin_1,zmin_2,tapperz, iflag)
!       end if
       
!       tapper = tapperx * tappery * tapperz

       
!       vp_final  = vp_final  * (1d0-tapper) + vp3d(ilon,ilat,idepth)  * (0.d0 +tapper)
!       vs_final  = vs_final  * (1.d0-tapper) + vs3d(ilon,ilat,idepth)  * (0.d0 +tapper)
!       rho_final = rho_final * (1.d0-tapper) + rho3d(ilon,ilat,idepth) * (0.d0 +tapper)
!!!!! end tapper

!    end if

  end subroutine model_3d_vm

!-----------------------------------------------
  subroutine tapper_3D( x, y, z, vp_final, vs_final, rho_final, vpnew, vsnew, rhonew)
    use external_model
    implicit none
    double precision :: tapperx, tappery, tapperz, tapper
    double precision :: vpnew,vsnew,rhonew
    double precision :: x, y, z
    real(kind=CUSTOM_REAL) rho_final,vp_final,vs_final
    integer izone,iflag

    call in_which_zone_I_am(x,y,z,izone)

    if (izone == 1) then
       vp_final  = vpnew       
       vs_final  = vsnew
       rho_final = rhonew
    end if

    if (izone == 2) then

       tapperx=1.d0
       tappery=1.d0
       tapperz=1.d0
       tapper=1.d0

       ! add tapper
       if     (x > xmin_1 .and. x < xmin_2) then
          iflag = 1
          call tapper_1D(x,xmin_1,xmin_2,tapperx, iflag)


       elseif (x < xmax_1 .and. x > xmax_2) then
          iflag = -1
          call tapper_1D(x,xmax_1,xmax_2,tapperx, iflag)

       end if

       if     (y > ymin_1 .and. y < ymin_2) then
          iflag = 1
          call tapper_1D(y,ymin_1,ymin_2,tappery, iflag)

       elseif (y < ymax_1 .and. y > ymax_2) then
          iflag = -1
          call tapper_1D(y,ymax_1,ymax_2,tappery, iflag)
       end if

       if     (z > zmin_1 .and. z < zmin_2) then
          iflag = 1
          call tapper_1D(z,zmin_1,zmin_2,tapperz, iflag)
       end if

       tapper = tapperx * tappery * tapperz

       vp_final  = vp_final  * (1.d0-tapper) + vpnew  * (0.d0 +tapper)
       vs_final  = vs_final  * (1.d0-tapper) + vsnew  * (0.d0 +tapper)
       rho_final = rho_final * (1.d0-tapper) + rhonew * (0.d0 +tapper)

   end if

  end subroutine tapper_3D


!----------------------------------------------
  subroutine tri_interp(i,j,k,ri, rj, rk , vpi,vsi, rhoi)
    use external_model
    implicit none
    integer i,j,k
    double precision ri,rj,rk,vpi,vsi, rhoi 
    double precision f,f1,f2,f3,f4,f5,f6,f7,f8
    double precision t,u,w

    t=ri
    u=rj
    w=rk

    !! vp
    f1=vp3d(i,j,k)
    f2=vp3d(i+1,j,k)
    f3=vp3d(i+1,j+1,k)
    f4=vp3d(i,j+1,k)
    f5=vp3d(i,j,k+1)
    f6=vp3d(i+1,j,k+1)
    f7=vp3d(i+1,j+1,k+1)
    f8=vp3d(i,j+1,k+1)
    f= (1.d0-t) * (1.d0-u) * (1.d0-w)*  f1+&
          t     * (1.d0-u) * (1.d0-w)*  f2+&
          t     *    u     * (1.d0-w)*  f3+&
       (1.d0-t) *    u     * (1.d0-w)*  f4+&
       (1.d0-t) * (1.d0-u) *   w     *  f5+&
          t     * (1.d0-u) *   w     *  f6+&
          t     *   u      *   w     *  f7+&
       (1.d0-t) *   u      *   w     *  f8

    vpi=f

    !! vs
    f1=vs3d(i,j,k)
    f2=vs3d(i+1,j,k)
    f3=vs3d(i+1,j+1,k)
    f4=vs3d(i,j+1,k)
    f5=vs3d(i,j,k+1)
    f6=vs3d(i+1,j,k+1)
    f7=vs3d(i+1,j+1,k+1)
    f8=vs3d(i,j+1,k+1)
    f= (1.d0-t) * (1.d0-u) * (1.d0-w)*  f1+&
          t     * (1.d0-u) * (1.d0-w)*  f2+&
          t     *    u     * (1.d0-w)*  f3+&
       (1.d0-t) *    u     * (1.d0-w)*  f4+&
       (1.d0-t) * (1.d0-u) *   w     *  f5+&
          t     * (1.d0-u) *   w     *  f6+&
          t     *   u      *   w     *  f7+&
       (1.d0-t) *   u      *   w     *  f8

    vsi=f

    !! rho
    f1=rho3d(i,j,k)
    f2=rho3d(i+1,j,k)
    f3=rho3d(i+1,j+1,k)
    f4=rho3d(i,j+1,k)
    f5=rho3d(i,j,k+1)
    f6=rho3d(i+1,j,k+1)
    f7=rho3d(i+1,j+1,k+1)
    f8=rho3d(i,j+1,k+1)
    f= (1.d0-t) * (1.d0-u) * (1.d0-w)*  f1+&
          t     * (1.d0-u) * (1.d0-w)*  f2+&
          t     *    u     * (1.d0-w)*  f3+&
       (1.d0-t) *    u     * (1.d0-w)*  f4+&
       (1.d0-t) * (1.d0-u) *   w     *  f5+&
          t     * (1.d0-u) *   w     *  f6+&
          t     *   u      *   w     *  f7+&
       (1.d0-t) *   u      *   w     *  f8

    rhoi=f


  end subroutine tri_interp

  subroutine model_2D(x_eval,y_eval,z_eval, &
                             rho_final,vp_final,vs_final,r1)
    use external_model
    implicit none
    double precision r1,radius
    double precision rho,vp,vs
    double precision x_eval,y_eval,z_eval
    double precision x,y,z
    real(kind=CUSTOM_REAL) rho_final,vp_final,vs_final
    double precision Interpol,Xtol
    integer idepth,ilon,izone,iflag
    double precision  :: tapperx, tappery, tapperz, tapper


!!$  !! HARDCODED ---------------
!!$    !! ZONE 1 externe
!!$    xmin_1=-185000.d0
!!$    xmax_1= 185000.d0
!!$
!!$    ymin_1=-79000.d0
!!$    ymax_1= 79000.d0
!!$    
!!$    zmin_1 = 5000.d0
!!$    
!!$    !! ZONE 2 intermediaire
!!$    xmin_2= -165000.d0
!!$    xmax_2=  165000.d0
!!$
!!$    ymin_2= -65000.d0
!!$    ymax_2=  65000.d0
!!$    
!!$    zmin_2 = 10000.d0

    x=x_eval
    y=y_eval
    z=z_eval
    Xtol=1d-2

    radius = dsqrt(x_eval**2 + y_eval**2 + (z_eval+zref)**2)
        
    ! 
    idepth = -floor((radius - (Radius_Earth - z2d(1,1) )) / h2d) + 1  
    ilon   = floor( (x_eval - x2d(1,1)) / h2d ) + 1

   
    if (idepth < 1) idepth = 1
    if (idepth > nz2d) idepth = nz2d

    if (ilon < 1) ilon = 1
    if (ilon > nx2d) ilon = nx2d
     
    
    
    call in_which_zone_I_am(x,y,z,izone) 


    if (izone == 1) then
       vp_final = vp2d(ilon,idepth)   ! * 1000.d0
       vs_final = vs2d(ilon,idepth)   ! * 1000.d0
       rho_final = rho2d(ilon,idepth) ! * 1000.d0
    end if
    
    if (izone == 2) then
       tapperx=1.d0
       tappery=1.d0
       tapperz=1.d0
       tapper=1.d0
       ! add tapper
       if     (x > xmin_1 .and. x < xmin_2) then
          iflag = 1
          call tapper_1D(x,xmin_1,xmin_2,tapperx, iflag)
          

       elseif (x < xmax_1 .and. x > xmax_2) then
          iflag = -1
          call tapper_1D(x,xmax_1,xmax_2,tapperx, iflag)

       end if
       
       if     (y > ymin_1 .and. y < ymin_2) then
          iflag = 1
          call tapper_1D(y,ymin_1,ymin_2,tappery, iflag)
          
       elseif (y < ymax_1 .and. y > ymax_2) then
          iflag = -1
          call tapper_1D(y,ymax_1,ymax_2,tappery, iflag)
       end if
       
       if     (z > zmin_1 .and. z < zmin_2) then
          iflag = 1
          call tapper_1D(z,zmin_1,zmin_2,tapperz, iflag)
       end if
       
       tapper = tapperx * tappery * tapperz
        !!to test 
       !tapper=1.d0
       vp_final  = vp_final  * (1d0-tapper) + vp2d(ilon,idepth)  * (0.d0 +tapper)
       vs_final  = vs_final  * (1.d0-tapper) + vs2d(ilon,idepth)  * (0.d0 +tapper)
       rho_final = rho_final * (1.d0-tapper) + rho2d(ilon,idepth) * (0.d0 +tapper)
    end if
    
    !write(*,*) idepth,radius,Radius_Earth-z2d(1,1),vp_final
  end subroutine model_2D



  subroutine in_which_zone_I_am(x,y,z,izone)
    use external_model
    implicit none
    double precision :: x, y, z
    integer izone

    izone = 3 ! par defaut on est dans la zone externe
 
    ! 
    if ( x > xmin_2 .and. x < xmax_2 .and. &
         y > ymin_2 .and. y < ymax_2 .and. & 
         z > zmin_2 ) then
       izone = 1 ! zone interne
       return
    end if
    !
     if ( x > xmin_1 .and. x < xmax_1 .and. &
         y > ymin_1 .and. y < ymax_1 .and. & 
         z > zmin_1 ) then
       izone = 2 ! zone intermediaire 
       return
    end if

  end subroutine in_which_zone_I_am

  subroutine tapper_1D(x,x0,x1,tapper, iflag)
    implicit none
    double precision, parameter ::  my_pi=3.141592653589793d0
    double precision :: te1, te2, tapper
    double precision :: x, x0, x1
    integer :: iflag

    if (iflag==-1) then
       te1 = my_pi/(x0 - x1)
       !te2=  -te1 * x0
       te2=  -te1 * x1
    elseif (iflag==1) then
       te1 = my_pi/(x1 - x0)
       !te2=  -te1 * x0
       te2=  -te1 * x1
    end if

    tapper = (0.5d0+0.5d0*cos(te1*x+te2))
    
  end subroutine tapper_1D
