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

  subroutine compute_arrays_source(ispec_selected_source, &
             xi_source,eta_source,gamma_source,sourcearray, &
             Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
             xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
             xigll,yigll,zigll,nspec)

  implicit none

  include "constants.h"

  integer ispec_selected_source
  integer nspec

  double precision xi_source,eta_source,gamma_source
  double precision Mxx,Myy,Mzz,Mxy,Mxz,Myz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xix,xiy,xiz,etax,etay,etaz, &
        gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLY,NGLLZ) :: sourcearray

  double precision xixd,xiyd,xizd,etaxd,etayd,etazd,gammaxd,gammayd,gammazd

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

! source arrays
  double precision, dimension(:,:,:,:), allocatable :: sourcearrayd
  double precision, dimension(:,:,:), allocatable :: G11,G12,G13,G21,G22,G23,G31,G32,G33
  double precision, dimension(:), allocatable :: hxis,hpxis,hetas,hpetas,hgammas,hpgammas

  integer k,l,m
  integer ir,it,iv

! allocate memory for arrays
  allocate(sourcearrayd(3,NGLLX,NGLLY,NGLLZ))
  allocate(G11(NGLLX,NGLLY,NGLLZ))
  allocate(G12(NGLLX,NGLLY,NGLLZ))
  allocate(G13(NGLLX,NGLLY,NGLLZ))
  allocate(G21(NGLLX,NGLLY,NGLLZ))
  allocate(G22(NGLLX,NGLLY,NGLLZ))
  allocate(G23(NGLLX,NGLLY,NGLLZ))
  allocate(G31(NGLLX,NGLLY,NGLLZ))
  allocate(G32(NGLLX,NGLLY,NGLLZ))
  allocate(G33(NGLLX,NGLLY,NGLLZ))
  allocate(hxis(NGLLX))
  allocate(hpxis(NGLLX))
  allocate(hetas(NGLLY))
  allocate(hpetas(NGLLY))
  allocate(hgammas(NGLLZ))
  allocate(hpgammas(NGLLZ))

! calculate G_ij for general source location
! the source does not necessarily correspond to a Gauss-Lobatto point
  do m=1,NGLLZ
    do l=1,NGLLY
      do k=1,NGLLX

        xixd    = dble(xix(k,l,m,ispec_selected_source))
        xiyd    = dble(xiy(k,l,m,ispec_selected_source))
        xizd    = dble(xiz(k,l,m,ispec_selected_source))
        etaxd   = dble(etax(k,l,m,ispec_selected_source))
        etayd   = dble(etay(k,l,m,ispec_selected_source))
        etazd   = dble(etaz(k,l,m,ispec_selected_source))
        gammaxd = dble(gammax(k,l,m,ispec_selected_source))
        gammayd = dble(gammay(k,l,m,ispec_selected_source))
        gammazd = dble(gammaz(k,l,m,ispec_selected_source))

        G11(k,l,m) = Mxx*xixd+Mxy*xiyd+Mxz*xizd
        G12(k,l,m) = Mxx*etaxd+Mxy*etayd+Mxz*etazd
        G13(k,l,m) = Mxx*gammaxd+Mxy*gammayd+Mxz*gammazd
        G21(k,l,m) = Mxy*xixd+Myy*xiyd+Myz*xizd
        G22(k,l,m) = Mxy*etaxd+Myy*etayd+Myz*etazd
        G23(k,l,m) = Mxy*gammaxd+Myy*gammayd+Myz*gammazd
        G31(k,l,m) = Mxz*xixd+Myz*xiyd+Mzz*xizd
        G32(k,l,m) = Mxz*etaxd+Myz*etayd+Mzz*etazd
        G33(k,l,m) = Mxz*gammaxd+Myz*gammayd+Mzz*gammazd

      enddo
    enddo
  enddo

! compute Lagrange polynomials at the source location
  call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  call lagrange_any(eta_source,NGLLY,yigll,hetas,hpetas)
  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)

! calculate source array
  do m=1,NGLLZ
    do l=1,NGLLY
      do k=1,NGLLX

      sourcearrayd(:,k,l,m) = ZERO

      do iv=1,NGLLZ
        do it=1,NGLLY
          do ir=1,NGLLX

        sourcearrayd(1,k,l,m) = sourcearrayd(1,k,l,m) + hxis(ir)*hetas(it)*hgammas(iv) &
                               *(G11(ir,it,iv)*hpxis(k)*hetas(l)*hgammas(m) &
                               +G12(ir,it,iv)*hxis(k)*hpetas(l)*hgammas(m) &
                               +G13(ir,it,iv)*hxis(k)*hetas(l)*hpgammas(m))

        sourcearrayd(2,k,l,m) = sourcearrayd(2,k,l,m) + hxis(ir)*hetas(it)*hgammas(iv) &
                               *(G21(ir,it,iv)*hpxis(k)*hetas(l)*hgammas(m) &
                               +G22(ir,it,iv)*hxis(k)*hpetas(l)*hgammas(m) &
                               +G23(ir,it,iv)*hxis(k)*hetas(l)*hpgammas(m))

        sourcearrayd(3,k,l,m) = sourcearrayd(3,k,l,m) + hxis(ir)*hetas(it)*hgammas(iv) &
                               *(G31(ir,it,iv)*hpxis(k)*hetas(l)*hgammas(m) &
                               +G32(ir,it,iv)*hxis(k)*hpetas(l)*hgammas(m) &
                               +G33(ir,it,iv)*hxis(k)*hetas(l)*hpgammas(m))

          enddo
        enddo
      enddo

      enddo
    enddo
  enddo

! distinguish whether single or double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then
    sourcearray(:,:,:,:) = sngl(sourcearrayd(:,:,:,:))
  else
    sourcearray(:,:,:,:) = sourcearrayd(:,:,:,:)
  endif

! deallocate arrays
  deallocate(sourcearrayd)
  deallocate(G11)
  deallocate(G12)
  deallocate(G13)
  deallocate(G21)
  deallocate(G22)
  deallocate(G23)
  deallocate(G31)
  deallocate(G32)
  deallocate(G33)
  deallocate(hxis)
  deallocate(hpxis)
  deallocate(hetas)
  deallocate(hpetas)
  deallocate(hgammas)
  deallocate(hpgammas)

  end subroutine compute_arrays_source

