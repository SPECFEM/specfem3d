!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 2
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine calc_jacobian(myrank,xixstore,xiystore,xizstore, &
     etaxstore,etaystore,etazstore, &
     gammaxstore,gammaystore,gammazstore,jacobianstore, &
     xstore,ystore,zstore,xelm,yelm,zelm,shape3D,dershape3D,ispec,nspec)

  implicit none

  include "constants.h"

  integer ispec,nspec,myrank

  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

  double precision, dimension(NGNOD) :: xelm,yelm,zelm

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore,jacobianstore

  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  integer i,j,k,ia
  double precision xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision xmesh,ymesh,zmesh
  double precision xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision jacobian

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

      xxi = ZERO
      xeta = ZERO
      xgamma = ZERO
      yxi = ZERO
      yeta = ZERO
      ygamma = ZERO
      zxi = ZERO
      zeta = ZERO
      zgamma = ZERO
      xmesh = ZERO
      ymesh = ZERO
      zmesh = ZERO

      do ia=1,NGNOD
        xxi = xxi + dershape3D(1,ia,i,j,k)*xelm(ia)
        xeta = xeta + dershape3D(2,ia,i,j,k)*xelm(ia)
        xgamma = xgamma + dershape3D(3,ia,i,j,k)*xelm(ia)
        yxi = yxi + dershape3D(1,ia,i,j,k)*yelm(ia)
        yeta = yeta + dershape3D(2,ia,i,j,k)*yelm(ia)
        ygamma = ygamma + dershape3D(3,ia,i,j,k)*yelm(ia)
        zxi = zxi + dershape3D(1,ia,i,j,k)*zelm(ia)
        zeta = zeta + dershape3D(2,ia,i,j,k)*zelm(ia)
        zgamma = zgamma + dershape3D(3,ia,i,j,k)*zelm(ia)
        xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
        ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
        zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
      enddo

      jacobian = xxi*(yeta*zgamma-ygamma*zeta) - &
             xeta*(yxi*zgamma-ygamma*zxi) + &
             xgamma*(yxi*zeta-yeta*zxi)

! can ignore negative jacobian in mesher if needed when debugging code
      if(jacobian <= ZERO) call exit_MPI(myrank,'3D Jacobian undefined')

!     invert the relation (Fletcher p. 50 vol. 2)
      xix = (yeta*zgamma-ygamma*zeta) / jacobian
      xiy = (xgamma*zeta-xeta*zgamma) / jacobian
      xiz = (xeta*ygamma-xgamma*yeta) / jacobian
      etax = (ygamma*zxi-yxi*zgamma) / jacobian
      etay = (xxi*zgamma-xgamma*zxi) / jacobian
      etaz = (xgamma*yxi-xxi*ygamma) / jacobian
      gammax = (yxi*zeta-yeta*zxi) / jacobian
      gammay = (xeta*zxi-xxi*zeta) / jacobian
      gammaz = (xxi*yeta-xeta*yxi) / jacobian

!     compute and store the jacobian for the solver
      jacobian = 1. / (xix*(etay*gammaz-etaz*gammay) &
                      -xiy*(etax*gammaz-etaz*gammax) &
                      +xiz*(etax*gammay-etay*gammax))

!     save the derivatives and the jacobian

! distinguish between single and double precision for reals
      if(CUSTOM_REAL == SIZE_REAL) then
        xixstore(i,j,k,ispec) = sngl(xix)
        xiystore(i,j,k,ispec) = sngl(xiy)
        xizstore(i,j,k,ispec) = sngl(xiz)
        etaxstore(i,j,k,ispec) = sngl(etax)
        etaystore(i,j,k,ispec) = sngl(etay)
        etazstore(i,j,k,ispec) = sngl(etaz)
        gammaxstore(i,j,k,ispec) = sngl(gammax)
        gammaystore(i,j,k,ispec) = sngl(gammay)
        gammazstore(i,j,k,ispec) = sngl(gammaz)
        jacobianstore(i,j,k,ispec) = sngl(jacobian)
      else
        xixstore(i,j,k,ispec) = xix
        xiystore(i,j,k,ispec) = xiy
        xizstore(i,j,k,ispec) = xiz
        etaxstore(i,j,k,ispec) = etax
        etaystore(i,j,k,ispec) = etay
        etazstore(i,j,k,ispec) = etaz
        gammaxstore(i,j,k,ispec) = gammax
        gammaystore(i,j,k,ispec) = gammay
        gammazstore(i,j,k,ispec) = gammaz
        jacobianstore(i,j,k,ispec) = jacobian
      endif

      xstore(i,j,k,ispec) = xmesh
      ystore(i,j,k,ispec) = ymesh
      zstore(i,j,k,ispec) = zmesh

      enddo
    enddo
  enddo

  end subroutine calc_jacobian

