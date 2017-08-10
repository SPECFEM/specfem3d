module interp_mod

  use global_parameters, only: CUSTOM_REAL,NGNOD,NGLLX,NGLLY
  real(kind=CUSTOM_REAL)    :: a,b,c,d,det_jacobian
  real(kind=CUSTOM_REAL) :: hxir(NGLLX),hetar(NGLLY)
  real(kind=CUSTOM_REAL) :: hpxir(NGLLX),hpetar(NGLLY)
  real(kind=CUSTOM_REAL) :: xigll(NGLLX),yigll(NGLLY)
  real(kind=CUSTOM_REAL) :: wxgll(NGLLX),wygll(NGLLY)
  real(kind=CUSTOM_REAL), parameter :: GAUSSALPHA=0._CUSTOM_REAL,GAUSSBETA=0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: zero=0._CUSTOM_REAL,one=1._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: three=3._CUSTOM_REAL,quart=0.25_CUSTOM_REAL,half=0.5_CUSTOM_REAL
contains

  subroutine interpole_field(xi,eta,field,interp_value)

    real(kind=CUSTOM_REAL) xi,eta,interp_value,hlagrange
    real(kind=CUSTOM_REAL) field(NGLLX,NGLLY)
    integer igll,jgll
    !! field(igll,jgll), nodes_gll(igll,jgll) !! pour un element

    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)

    call lagrange_any(xi,NGLLX,xigll,hxir,hpxir)
    call lagrange_any(eta,NGLLY,yigll,hetar,hpetar)


    interp_value=0.d0


    do jgll = 1,NGLLY
       do igll = 1,NGLLX

          ! ploynome de lagrange
          hlagrange = hxir(igll)*hetar(jgll) !*hgammar(kgll)

          ! interpolation de la valeur
          interp_value = interp_value + field(igll,jgll)*hlagrange

       enddo
    enddo

  end subroutine interpole_field

  subroutine find_xix_eta(nodes_crd,xi,eta,s_target,z_target)

    real(kind=CUSTOM_REAL)    :: xi,eta,s,z,s_target,z_target
    real(kind=CUSTOM_REAL)    :: nodes_crd(ngnod,2)

    real(kind=CUSTOM_REAL)    :: sph(ngnod)
    real(kind=CUSTOM_REAL)    :: distmin,dist
    real(kind=CUSTOM_REAL)    :: dxi,deta
    integer inode,iguess,iter_newton,niter_newton

    !zero=0.;one=1.;quart=0.25;half=0.5;three=3.
    niter_newton=6
    ! find the closest node
    distmin=1d30
    iguess=1
    do inode = 1, ngnod
       dist=(nodes_crd(inode,1)-s_target)**2 + (nodes_crd(inode,2)-z_target)**2
       if (dist < distmin) then
          distmin=dist
          iguess=inode
       endif
    enddo

    ! convert to xi,eta initial guess  !! base sur le dessin dessous (a verifier)
    if (iguess == 1) then
       xi=-1.
       eta=-1.
    endif
    if (iguess == 2) then
       xi=0.
       eta=-1.
    endif
    if (iguess == 3) then
       xi=1.
       eta=-1.
    endif
    if (iguess == 4) then
       xi=1.
       eta=0.
    endif
    if (iguess == 5) then
       xi=1.
       eta=1.
    endif
    if (iguess == 6) then
       xi=0.
       eta=1.
    endif
    if (iguess == 7) then
       xi=-1.
       eta=1.
    endif
    if (iguess == 8) then
       xi=-1.
       eta=0.
    endif
    !write(*,*) xi,eta
    do iter_newton=1,niter_newton
       det_jacobian =  det_jacobian_shape(xi, eta, nodes_crd)
       call shp8(xi,eta,sph)
       s=0.
       z=0.
       do inode=1,ngnod
          s=s+sph(inode)*nodes_crd(inode,1)
          z=z+sph(inode)*nodes_crd(inode,2)
       enddo

       dxi=(d*(s-s_target)  -b  *(z-z_target))/det_jacobian
       deta=(-c*(s-s_target) +a  *(z-z_target))/det_jacobian

       xi  = xi  - dxi
       eta = eta - deta
       !write(*,*) xi,eta
    enddo
    !call det_jacobian_shape(xi, eta, nodes_crd,a,b,c,d)
    call shp8(xi,eta,sph)
    s=0.
    z=0.
    do inode=1,ngnod
       s=s+sph(inode)*nodes_crd(inode,1)
       z=z+sph(inode)*nodes_crd(inode,2)
    enddo
!!$    write(*,*) 'distance :',s-s_target,z-z_target
!!$    write(*,*) s,s_target
!!$    write(*,*) z,z_target
!!$    write(*,*) 'xi, eta :', xi,eta
!!$    write(*,*)

  end subroutine find_xix_eta

  !-----------------------------------------------------------------------------------------
  subroutine shp8(xil,etal,shp)
    !
    ! This routine computes and returns the quadratic
    ! shape functions axixiociated with a 8-nodes serendip
    ! element for a given point of coordinates (xi,eta).
    !
    ! Topology is defined as follows
    !
    ! 7 - - - 6 - - - 5
    ! |       ^       |
    ! |   eta |       |
    ! |       |       |
    ! 8        --->   4
    ! |        xi     |
    ! |               |
    ! |               |
    ! 1 - - - 2 - - - 3
    !

    real(kind=CUSTOM_REAL)    :: xil, etal
    real(kind=CUSTOM_REAL)    :: shp(8)
    real(kind=CUSTOM_REAL)    :: xip,xim,etap,etam,xixi,etaeta

    shp(:) = zero


    xip    = one +  xil
    xim    = one -  xil
    etap   = one + etal
    etam   = one - etal
    xixi   =  xil *  xil
    etaeta = etal * etal

    ! Corners first:
    shp(1) = quart * xim * etam * (xim + etam - three)
    shp(3) = quart * xip * etam * (xip + etam - three)
    shp(5) = quart * xip * etap * (xip + etap - three)
    shp(7) = quart * xim * etap * (xim + etap - three)

    ! Then midpoints:
    shp(2) = half  * etam * (one -   xixi)
    shp(4) = half  *  xip * (one - etaeta)
    shp(6) = half  * etap * (one -   xixi)
    shp(8) = half  *  xim * (one - etaeta)

  end subroutine shp8
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  subroutine shp8der(xil,etal,shpder)
    !
    ! This routine computes and returns the derivatives
    ! of the shape functions axixiociated with a 8-nodes serendip
    ! element for a given point of coordinates (xi,eta).
    !
    ! Topology is defined as follows
    !
    ! 7 - - - 6 - - - 5
    ! |       ^       |
    ! |   eta |       |
    ! |       |       |
    ! 8        --->   4
    ! |        xi     |
    ! |               |
    ! |               |
    ! 1 - - - 2 - - - 3
    !
    !
    ! shpder(:,1) : derivative wrt xi
    ! shpder(:,2) : derivative wrt eta

    real(kind=CUSTOM_REAL)    :: xil, etal
    real(kind=CUSTOM_REAL)    :: shpder(8,2)
    real(kind=CUSTOM_REAL)    :: xip,xim,etap,etam,xixi,etaeta

    shpder(:,:) = zero

    xip    = one +  xil
    xim    = one -  xil
    etap   = one + etal
    etam   = one - etal
    xixi   =  xil *  xil
    etaeta = etal * etal

  ! Corners first:
    shpder(1,1) = -quart * etam * ( xim + xim + etam - three)
    shpder(1,2) = -quart *  xim * (etam + xim + etam - three)
    shpder(3,1) =  quart * etam * ( xip + xip + etam - three)
    shpder(3,2) = -quart *  xip * (etam + xip + etam - three)
    shpder(5,1) =  quart * etap * ( xip + xip + etap - three)
    shpder(5,2) =  quart *  xip * (etap + xip + etap - three)
    shpder(7,1) = -quart * etap * ( xim + xim + etap - three)
    shpder(7,2) =  quart *  xim * (etap + xim + etap - three)

    ! Then midside points :
    shpder(2,1) = -one  * xil * etam
    shpder(2,2) = -half * (one - xixi)
    shpder(4,1) =  half * (one - etaeta)
    shpder(4,2) = -one  * etal * xip
    shpder(6,1) = -one  * xil * etap
    shpder(6,2) =  half * (one - xixi)
    shpder(8,1) = -half * (one - etaeta)
    shpder(8,2) = -one  * etal * xim

  end subroutine shp8der
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  real(kind=CUSTOM_REAL)    function det_jacobian_shape(xil, etal, nodes_crd)
    ! This routines the value of the Jacobian (that is,
    ! the determinant of the Jacobian matrix), for any point
    ! inside a given element. IT ASSUMES 8 nodes 2D isoparametric
    ! formulation of the geometrical transformation and therefore
    ! requires the knowledge of the coordinated of the 8 control
    ! points, which are defined as follows :
    !
    !     7 - - - 6 - - - 5
    !     |       ^       |
    !     |   eta |       |
    !     |       |       |
    !     8        --->   4
    !     |        xi     |
    !     |               |
    !     |               |
    !     1 - - - 2 - - - 3 .
    !
    real(kind=CUSTOM_REAL)    :: xil, etal, nodes_crd(8,2)
    integer :: inode
    real(kind=CUSTOM_REAL)    :: shpder(8,2)!,a,d,b,c

    ! Compute the appropriate derivatives of the shape
    ! functions

    call shp8der(xil,etal,shpder)

    a = zero
    b = zero
    c = zero
    d = zero

    do inode = 1, 8
       a = a + nodes_crd(inode,1)*shpder(inode,1)
       d = d + nodes_crd(inode,2)*shpder(inode,2)
       b = b + nodes_crd(inode,1)*shpder(inode,2)
       c = c + nodes_crd(inode,2)*shpder(inode,1)
    enddo

    det_jacobian_shape = a*d - b*c


  end function det_jacobian_shape


  subroutine lagrange_any(xi,NGLL,xigll,h,hprime)

! subroutine to compute the Lagrange interpolants based upon the GLL points
! and their first derivatives at any point xi in [-1,1]

  implicit none

  integer, intent(in) :: NGLL
  real(kind=CUSTOM_REAL), intent(in) :: xi,xigll(NGLL)
  real(kind=CUSTOM_REAL), intent(out) :: h(NGLL),hprime(NGLL)

  integer dgr,i,j
  real(CUSTOM_REAL) prod1,prod2

  do dgr=1,NGLL

    prod1 = 1.0d0
    prod2 = 1.0d0
    do i=1,NGLL
      if (i /= dgr) then
        prod1 = prod1*(xi-xigll(i))
        prod2 = prod2*(xigll(dgr)-xigll(i))
      endif
    enddo
    h(dgr)=prod1/prod2

    hprime(dgr)=0.0d0
    do i=1,NGLL
      if (i /= dgr) then
        prod1=1.0d0
        do j=1,NGLL
          if (j /= dgr .and. j /= i) prod1 = prod1*(xi-xigll(j))
        enddo
        hprime(dgr) = hprime(dgr)+prod1
      endif
    enddo
    hprime(dgr) = hprime(dgr)/prod2

  enddo

  end subroutine lagrange_any

end module interp_mod
