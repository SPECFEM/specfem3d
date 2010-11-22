!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
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

  subroutine calc_jacobian(myrank,xixstore,xiystore,xizstore, &
                          etaxstore,etaystore,etazstore, &
                          gammaxstore,gammaystore,gammazstore,jacobianstore, &
                          xstore,ystore,zstore, &
                          xelm,yelm,zelm,shape3D,dershape3D,ispec,nspec)

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

!
!-------------------------------------------------------------------------------------------------
!

! This subroutine recomputes the 3D jacobian for one element
! based upon all GLL points
! Hejun Zhu OCT16,2009

! input: myrank,
!        xstore,ystore,zstore ----- input position
!        xigll,yigll,zigll ----- gll points position
!        ispec,nspec       ----- element number
!        ACTUALLY_STORE_ARRAYS   ------ save array or not

! output: xixstore,xiystore,xizstore,
!         etaxstore,etaystore,etazstore,
!         gammaxstore,gammaystore,gammazstore ------ parameters used for calculating jacobian
!
!
!  subroutine recalc_jacobian_gll3D(myrank,xixstore,xiystore,xizstore, &
!                                  etaxstore,etaystore,etazstore, &
!                                  gammaxstore,gammaystore,gammazstore,jacobianstore, &
!                                  xstore,ystore,zstore, &
!                                  ispec,nspec, &
!                                  xigll,yigll,zigll, &
!                                  ACTUALLY_STORE_ARRAYS)
!
!  implicit none
!
!  include "constants.h"
!
!  ! input parameter
!  integer::myrank,ispec,nspec
!  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore
!  double precision, dimension(NGLLX):: xigll
!  double precision, dimension(NGLLY):: yigll
!  double precision, dimension(NGLLZ):: zigll
!  logical::ACTUALLY_STORE_ARRAYS
!
!
!  ! output results
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
!                        xixstore,xiystore,xizstore,&
!                        etaxstore,etaystore,etazstore,&
!                        gammaxstore,gammaystore,gammazstore,&
!                        jacobianstore
!
!
!  ! other parameters for this subroutine
!  integer:: i,j,k,i1,j1,k1
!  double precision:: xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
!  double precision:: xi,eta,gamma
!  double precision,dimension(NGLLX):: hxir,hpxir
!  double precision,dimension(NGLLY):: hetar,hpetar
!  double precision,dimension(NGLLZ):: hgammar,hpgammar
!  double precision:: hlagrange,hlagrange_xi,hlagrange_eta,hlagrange_gamma
!  double precision:: jacobian
!  double precision:: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
!
!
!
!  ! test parameters which can be deleted
!  double precision:: xmesh,ymesh,zmesh
!  double precision:: sumshape,sumdershapexi,sumdershapeeta,sumdershapegamma
!
!  ! first go over all 125 gll points
!  do k=1,NGLLZ
!    do j=1,NGLLY
!      do i=1,NGLLX
!
!            xxi = 0.0
!            xeta = 0.0
!            xgamma = 0.0
!            yxi = 0.0
!            yeta = 0.0
!            ygamma = 0.0
!            zxi = 0.0
!            zeta = 0.0
!            zgamma = 0.0
!
!            xi = xigll(i)
!            eta = yigll(j)
!            gamma = zigll(k)
!
!            ! calculate lagrange polynomial and its derivative
!            call lagrange_any(xi,NGLLX,xigll,hxir,hpxir)
!            call lagrange_any(eta,NGLLY,yigll,hetar,hpetar)
!            call lagrange_any(gamma,NGLLZ,zigll,hgammar,hpgammar)
!
!            ! test parameters
!            sumshape = 0.0
!            sumdershapexi = 0.0
!            sumdershapeeta = 0.0
!            sumdershapegamma = 0.0
!            xmesh = 0.0
!            ymesh = 0.0
!            zmesh = 0.0
!
!
!            do k1 = 1,NGLLZ
!               do j1 = 1,NGLLY
!                  do i1 = 1,NGLLX
!                     hlagrange = hxir(i1)*hetar(j1)*hgammar(k1)
!                     hlagrange_xi = hpxir(i1)*hetar(j1)*hgammar(k1)
!                     hlagrange_eta = hxir(i1)*hpetar(j1)*hgammar(k1)
!                     hlagrange_gamma = hxir(i1)*hetar(j1)*hpgammar(k1)
!
!
!                     xxi = xxi + xstore(i1,j1,k1,ispec)*hlagrange_xi
!                     xeta = xeta + xstore(i1,j1,k1,ispec)*hlagrange_eta
!                     xgamma = xgamma + xstore(i1,j1,k1,ispec)*hlagrange_gamma
!
!                     yxi = yxi + ystore(i1,j1,k1,ispec)*hlagrange_xi
!                     yeta = yeta + ystore(i1,j1,k1,ispec)*hlagrange_eta
!                     ygamma = ygamma + ystore(i1,j1,k1,ispec)*hlagrange_gamma
!
!                     zxi = zxi + zstore(i1,j1,k1,ispec)*hlagrange_xi
!                     zeta = zeta + zstore(i1,j1,k1,ispec)*hlagrange_eta
!                     zgamma = zgamma + zstore(i1,j1,k1,ispec)*hlagrange_gamma
!
!                     ! test the lagrange polynomial and its derivate
!                     xmesh = xmesh + xstore(i1,j1,k1,ispec)*hlagrange
!                     ymesh = ymesh + ystore(i1,j1,k1,ispec)*hlagrange
!                     zmesh = zmesh + zstore(i1,j1,k1,ispec)*hlagrange
!                     sumshape = sumshape + hlagrange
!                     sumdershapexi = sumdershapexi + hlagrange_xi
!                     sumdershapeeta = sumdershapeeta + hlagrange_eta
!                     sumdershapegamma = sumdershapegamma + hlagrange_gamma
!
!                  end do
!               end do
!            end do
!
!            ! Check the lagrange polynomial and its derivative
!            if (xmesh /=xstore(i,j,k,ispec).or.ymesh/=ystore(i,j,k,ispec).or.zmesh/=zstore(i,j,k,ispec)) then
!                    call exit_MPI(myrank,'new mesh positions are wrong in recalc_jacobian_gall3D.f90')
!            end if
!            if(abs(sumshape-one) >  TINYVAL) then
!                    call exit_MPI(myrank,'error shape functions in recalc_jacobian_gll3D.f90')
!            end if
!            if(abs(sumdershapexi) >  TINYVAL) then
!                    call exit_MPI(myrank,'error derivative xi shape functions in recalc_jacobian_gll3D.f90')
!            end if
!            if(abs(sumdershapeeta) >  TINYVAL) then
!                    call exit_MPI(myrank,'error derivative eta shape functions in recalc_jacobian_gll3D.f90')
!            end if
!            if(abs(sumdershapegamma) >  TINYVAL) then
!                    call exit_MPI(myrank,'error derivative gamma shape functions in recalc_jacobian_gll3D.f90')
!            end if
!
!
!            jacobian = xxi*(yeta*zgamma-ygamma*zeta) - &
!                 xeta*(yxi*zgamma-ygamma*zxi) + &
!                 xgamma*(yxi*zeta-yeta*zxi)
!
!            ! Check the jacobian
!            if(jacobian <= ZERO) then
!                   call exit_MPI(myrank,'3D Jacobian undefined in recalc_jacobian_gll3D.f90')
!            end if
!
!            !     invert the relation (Fletcher p. 50 vol. 2)
!            xix = (yeta*zgamma-ygamma*zeta) / jacobian
!            xiy = (xgamma*zeta-xeta*zgamma) / jacobian
!            xiz = (xeta*ygamma-xgamma*yeta) / jacobian
!            etax = (ygamma*zxi-yxi*zgamma) / jacobian
!            etay = (xxi*zgamma-xgamma*zxi) / jacobian
!            etaz = (xgamma*yxi-xxi*ygamma) / jacobian
!            gammax = (yxi*zeta-yeta*zxi) / jacobian
!            gammay = (xeta*zxi-xxi*zeta) / jacobian
!            gammaz = (xxi*yeta-xeta*yxi) / jacobian
!
!
!            !     compute and store the jacobian for the solver
!            jacobian = 1. / (xix*(etay*gammaz-etaz*gammay) &
!                            -xiy*(etax*gammaz-etaz*gammax) &
!                            +xiz*(etax*gammay-etay*gammax))
!
!            ! resave the derivatives and the jacobian
!            ! distinguish between single and double precision for reals
!            if (ACTUALLY_STORE_ARRAYS) then
!
!                if (myrank == 0) then
!                        print*,'xix before',xixstore(i,j,k,ispec),'after',xix
!                        print*,'etax before',etaxstore(i,j,k,ispec),'after',etax
!                        print*,'gammax before',gammaxstore(i,j,k,ispec),'after',gammax
!                end if
!
!                if(CUSTOM_REAL == SIZE_REAL) then
!                    xixstore(i,j,k,ispec) = sngl(xix)
!                    xiystore(i,j,k,ispec) = sngl(xiy)
!                    xizstore(i,j,k,ispec) = sngl(xiz)
!                    etaxstore(i,j,k,ispec) = sngl(etax)
!                    etaystore(i,j,k,ispec) = sngl(etay)
!                    etazstore(i,j,k,ispec) = sngl(etaz)
!                    gammaxstore(i,j,k,ispec) = sngl(gammax)
!                    gammaystore(i,j,k,ispec) = sngl(gammay)
!                    gammazstore(i,j,k,ispec) = sngl(gammaz)
!                    jacobianstore(i,j,k,ispec) = sngl(jacobian)
!                else
!                    xixstore(i,j,k,ispec) = xix
!                    xiystore(i,j,k,ispec) = xiy
!                    xizstore(i,j,k,ispec) = xiz
!                    etaxstore(i,j,k,ispec) = etax
!                    etaystore(i,j,k,ispec) = etay
!                    etazstore(i,j,k,ispec) = etaz
!                    gammaxstore(i,j,k,ispec) = gammax
!                    gammaystore(i,j,k,ispec) = gammay
!                    gammazstore(i,j,k,ispec) = gammaz
!                    jacobianstore(i,j,k,ispec) = jacobian
!                endif
!             end if
!        enddo
!    enddo
!  enddo
!
!  end subroutine recalc_jacobian_gll3D
!
