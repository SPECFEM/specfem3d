! This module outputs gravity field at required locations
!
! Authors:
! Surendra Somala (Caltech) surendra@caltech.edu - 2013
! with advice from Jan Harms and Pablo Ampuero (ampuero@gps.caltech.edu)

module gravity_perturbation

  implicit none

  include 'constants.h'

  private

  integer nstat,ntimgap,nstat_local
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstat,ystat,zstat
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accE,accN,accZ
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rho0_wm

  logical, save :: GRAVITY_SIMULATION = .false.

  public :: gravity_init, gravity_timeseries, gravity_output, GRAVITY_SIMULATION

contains

!=====================================================================

subroutine gravity_init()

  use specfem_par, only : NGLOB_AB, NSTEP, NSPEC_AB, mustore, &
       xstore, ystore, zstore, &
       xigll, yigll, zigll, &
       wxgll, wygll, wzgll, &
       NGNOD, ibool, myrank, IMAIN
  use specfem_par_elastic, only : rho_vs
  implicit none

  integer, parameter :: IIN_G = 367
  integer ier
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ):: rho_elem,vs_elem
  integer :: i,j,k,iglob,ispec,istat
  double precision :: Jac3D
  ! coordinates of the control points
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)
  integer ia
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz,iax,iay,iaz
  integer nstep_grav

  open(unit=IIN_G,file='../DATA/gravity_stations',status='old',iostat=ier)
  if( ier /= 0 ) then
    ! user output
    if( myrank == 0 ) then
      write(IMAIN,*)
      write(IMAIN,*) 'no gravity simulation'
      write(IMAIN,*)
    endif

    return
  endif

  GRAVITY_SIMULATION = .true.

  read(IIN_G,*) nstat,ntimgap

  ! user output
  if( myrank == 0 ) then
    write(IMAIN,*)
    write(IMAIN,*) 'incorporating gravity simulation'
    write(IMAIN,*) '    gravity stations: ',nstat
    write(IMAIN,*)
  endif

  allocate(xstat(nstat))
  allocate(ystat(nstat))
  allocate(zstat(nstat))
  do istat=1,nstat
     read(IIN_G,*) xstat(istat),ystat(istat),zstat(istat)
  enddo
  close(IIN_G)

  nstep_grav = floor(dble(NSTEP/ntimgap))
  allocate(accE(nstep_grav,nstat))
  allocate(accN(nstep_grav,nstat))
  allocate(accZ(nstep_grav,nstat))

  accE = 0._CUSTOM_REAL
  accN = 0._CUSTOM_REAL
  accZ = 0._CUSTOM_REAL

  allocate(rho0_wm(NGLOB_AB))
  rho0_wm = 0._CUSTOM_REAL

  call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)

     ! define coordinates of the control points of the element
  do ia=1,NGNOD

     if(iaddx(ia) == 0) then
        iax(ia) = 1
     else if(iaddx(ia) == 1) then
        iax(ia) = (NGLLX+1)/2
     else if(iaddx(ia) == 2) then
        iax(ia) = NGLLX
     else
        call exit_MPI(myrank,'incorrect value of iaddx')
     endif

     if(iaddy(ia) == 0) then
        iay(ia) = 1
     else if(iaddy(ia) == 1) then
        iay(ia) = (NGLLY+1)/2
     else if(iaddy(ia) == 2) then
        iay(ia) = NGLLY
     else
        call exit_MPI(myrank,'incorrect value of iaddy')
     endif

     if(iaddz(ia) == 0) then
        iaz(ia) = 1
     else if(iaddz(ia) == 1) then
        iaz(ia) = (NGLLZ+1)/2
     else if(iaddz(ia) == 2) then
        iaz(ia) = NGLLZ
     else
        call exit_MPI(myrank,'incorrect value of iaddz')
     endif

  enddo

  do ispec=1,NSPEC_AB

     rho_elem = rho_vs(:,:,:,ispec)*rho_vs(:,:,:,ispec)/mustore(:,:,:,ispec)

     do ia=1,NGNOD
        iglob = ibool(iax(ia),iay(ia),iaz(ia),ispec)
        xelm(ia) = dble(xstore(iglob))
        yelm(ia) = dble(ystore(iglob))
        zelm(ia) = dble(zstore(iglob))
     enddo

     do k = 1,NGLLZ
        do j = 1,NGLLY
           do i = 1,NGLLX
              iglob = ibool(i,j,k,ispec)
              call recompute_jacobian_gravity(xelm,yelm,zelm,xigll(i),yigll(j),zigll(k),Jac3D)
              rho0_wm(iglob) = rho0_wm(iglob) + Jac3D * wxgll(i) * wygll(j) * wzgll(k) * rho_elem(i,j,k)
           enddo
        enddo
     enddo

  enddo

end subroutine gravity_init

!=====================================================================

! recompute 3D jacobian at a given point for a 8-node element : modified from recompute_jacobian

subroutine recompute_jacobian_gravity(xelm,yelm,zelm,xi,eta,gamma,jacobian)

  use specfem_par, only : NGNOD
  implicit none

  include "constants.h"

  double precision x,y,z
  double precision xi,eta,gamma,jacobian

!  integer NGNOD

! coordinates of the control points
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

! 3D shape functions and their derivatives at receiver
  double precision shape3D(NGNOD)
  double precision dershape3D(NDIM,NGNOD)

  double precision xxi,yxi,zxi
  double precision xeta,yeta,zeta
  double precision xgamma,ygamma,zgamma
  double precision ra1,ra2,rb1,rb2,rc1,rc2

  integer ia

! for 8-node element
  double precision, parameter :: ONE_EIGHTH = 0.125d0

! recompute jacobian for any (xi,eta,gamma) point, not necessarily a GLL point

! check that the parameter file is correct
  if(NGNOD /= 8 ) &
       stop 'elements should have 8  control nodes'
!  if(NGNOD == 8) then

! ***
! *** create the 3D shape functions and the Jacobian for an 8-node element
! ***

!--- case of an 8-node 3D element (Dhatt-Touzot p. 115)

    ra1 = one + xi
    ra2 = one - xi

    rb1 = one + eta
    rb2 = one - eta

    rc1 = one + gamma
    rc2 = one - gamma

    shape3D(1) = ONE_EIGHTH*ra2*rb2*rc2
    shape3D(2) = ONE_EIGHTH*ra1*rb2*rc2
    shape3D(3) = ONE_EIGHTH*ra1*rb1*rc2
    shape3D(4) = ONE_EIGHTH*ra2*rb1*rc2
    shape3D(5) = ONE_EIGHTH*ra2*rb2*rc1
    shape3D(6) = ONE_EIGHTH*ra1*rb2*rc1
    shape3D(7) = ONE_EIGHTH*ra1*rb1*rc1
    shape3D(8) = ONE_EIGHTH*ra2*rb1*rc1

    dershape3D(1,1) = - ONE_EIGHTH*rb2*rc2
    dershape3D(1,2) = ONE_EIGHTH*rb2*rc2
    dershape3D(1,3) = ONE_EIGHTH*rb1*rc2
    dershape3D(1,4) = - ONE_EIGHTH*rb1*rc2
    dershape3D(1,5) = - ONE_EIGHTH*rb2*rc1
    dershape3D(1,6) = ONE_EIGHTH*rb2*rc1
    dershape3D(1,7) = ONE_EIGHTH*rb1*rc1
    dershape3D(1,8) = - ONE_EIGHTH*rb1*rc1

    dershape3D(2,1) = - ONE_EIGHTH*ra2*rc2
    dershape3D(2,2) = - ONE_EIGHTH*ra1*rc2
    dershape3D(2,3) = ONE_EIGHTH*ra1*rc2
    dershape3D(2,4) = ONE_EIGHTH*ra2*rc2
    dershape3D(2,5) = - ONE_EIGHTH*ra2*rc1
    dershape3D(2,6) = - ONE_EIGHTH*ra1*rc1
    dershape3D(2,7) = ONE_EIGHTH*ra1*rc1
    dershape3D(2,8) = ONE_EIGHTH*ra2*rc1

    dershape3D(3,1) = - ONE_EIGHTH*ra2*rb2
    dershape3D(3,2) = - ONE_EIGHTH*ra1*rb2
    dershape3D(3,3) = - ONE_EIGHTH*ra1*rb1
    dershape3D(3,4) = - ONE_EIGHTH*ra2*rb1
    dershape3D(3,5) = ONE_EIGHTH*ra2*rb2
    dershape3D(3,6) = ONE_EIGHTH*ra1*rb2
    dershape3D(3,7) = ONE_EIGHTH*ra1*rb1
    dershape3D(3,8) = ONE_EIGHTH*ra2*rb1


! compute coordinates and jacobian matrix
  x=ZERO
  y=ZERO
  z=ZERO
  xxi=ZERO
  xeta=ZERO
  xgamma=ZERO
  yxi=ZERO
  yeta=ZERO
  ygamma=ZERO
  zxi=ZERO
  zeta=ZERO
  zgamma=ZERO

  do ia=1,NGNOD
    x=x+shape3D(ia)*xelm(ia)
    y=y+shape3D(ia)*yelm(ia)
    z=z+shape3D(ia)*zelm(ia)

    xxi=xxi+dershape3D(1,ia)*xelm(ia)
    xeta=xeta+dershape3D(2,ia)*xelm(ia)
    xgamma=xgamma+dershape3D(3,ia)*xelm(ia)
    yxi=yxi+dershape3D(1,ia)*yelm(ia)
    yeta=yeta+dershape3D(2,ia)*yelm(ia)
    ygamma=ygamma+dershape3D(3,ia)*yelm(ia)
    zxi=zxi+dershape3D(1,ia)*zelm(ia)
    zeta=zeta+dershape3D(2,ia)*zelm(ia)
    zgamma=zgamma+dershape3D(3,ia)*zelm(ia)
  enddo

  jacobian = xxi*(yeta*zgamma-ygamma*zeta) - xeta*(yxi*zgamma-ygamma*zxi) + xgamma*(yxi*zeta-yeta*zxi)

  if(jacobian <= ZERO) stop '3D Jacobian undefined'


end subroutine recompute_jacobian_gravity

!=====================================================================

subroutine gravity_timeseries()

  use specfem_par, only : xstore, ystore, zstore, it, NGLOB_AB
  use specfem_par_elastic, only : displ
  implicit none

  real(kind=CUSTOM_REAL) :: G_const = 6.674e-11_CUSTOM_REAL
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: accEdV,accNdV,accZdV
  real(kind=CUSTOM_REAL) :: E_local,N_local,Z_local,E_all,N_all,Z_all
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: Rg,dotP
  integer :: istat, it_grav

  if( mod(it,ntimgap)==0 ) then
     it_grav = nint(dble(it/ntimgap))
     allocate(Rg(NGLOB_AB))
     allocate(dotP(NGLOB_AB))

     do istat=1,nstat
        Rg = sqrt((xstore-xstat(istat))**2+(ystore-ystat(istat))**2+(zstore-zstat(istat))**2)
        dotP = (xstore-xstat(istat))*displ(1,:)+(ystore-ystat(istat))*displ(2,:)+(zstore-zstat(istat))*displ(3,:)

        accEdV = G_const*rho0_wm/Rg**3*(displ(1,:)-3._CUSTOM_REAL*(dotP*(xstore-xstat(istat)))/Rg**2)
        E_local = sum(accEdV(:))
        call sum_all_all_cr(E_local,E_all)
        accE(it_grav,istat) = E_all

        accNdV = G_const*rho0_wm/Rg**3*(displ(2,:)-3._CUSTOM_REAL*(dotP*(ystore-ystat(istat)))/Rg**2)
        N_local = sum(accNdV(:))
        call sum_all_all_cr(N_local,N_all)
        accN(it_grav,istat) = N_all

        accZdV = G_const*rho0_wm/Rg**3*(displ(3,:)-3._CUSTOM_REAL*(dotP*(zstore-zstat(istat)))/Rg**2)
        Z_local = sum(accZdV(:))
        call sum_all_all_cr(Z_local,Z_all)
        accZ(it_grav,istat) = Z_all
     enddo
  endif

end subroutine gravity_timeseries

!=====================================================================

subroutine gravity_output()

  use specfem_par, only : myrank,NPROC,NSTEP,DT
  implicit none

  integer :: isample,istat,nstep_grav
  character(len=150) sisname

  nstep_grav = floor(dble(NSTEP/ntimgap))
  nstat_local = nint(dble(nstat/NPROC))

  do istat=1,nstat
     if(istat < myrank*nstat_local+1 .or. istat > (myrank+1)*nstat_local) cycle
     write(sisname,"('../OUTPUT_FILES/stat',I0,'.grav')") istat
     open(unit=IOUT,file=sisname,status='replace')
     do isample = 1,nstep_grav
        write(IOUT,*) isample*DT*ntimgap, accE(isample,istat),accN(isample,istat),accZ(isample,istat)
     enddo
     close(IOUT)
  enddo

  if(myrank==0) then !left-over stations
     do istat=NPROC*nstat_local,nstat
        write(sisname,"('../OUTPUT_FILES/stat',I0,'.grav')") istat
        open(unit=IOUT,file=sisname,status='replace')
        do isample = 1,nstep_grav
           write(IOUT,*) isample*DT*ntimgap, accE(isample,istat),accN(isample,istat),accZ(isample,istat)
        enddo
        close(IOUT)
     enddo
  endif

end subroutine gravity_output

!=====================================================================


end module gravity_perturbation
