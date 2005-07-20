!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2005
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

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray

  double precision xixd,xiyd,xizd,etaxd,etayd,etazd,gammaxd,gammayd,gammazd

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: G11,G12,G13,G21,G22,G23,G31,G32,G33
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer k,l,m
  integer ir,it,iv

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

! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then
    sourcearray(:,:,:,:) = sngl(sourcearrayd(:,:,:,:))
  else
    sourcearray(:,:,:,:) = sourcearrayd(:,:,:,:)
  endif

  end subroutine compute_arrays_source


!=================================================

subroutine compute_arrays_adjoint_source(myrank, adj_source_file, &
      xi_receiver,eta_receiver,gamma_receiver, adj_sourcearray, &
      xigll,yigll,zigll,NSTEP)


  implicit none

  include 'constants.h'

! input
  integer myrank, NSTEP

  double precision xi_receiver, eta_receiver, gamma_receiver

  character(len=*) adj_source_file

! output
  real(kind=CUSTOM_REAL) :: adj_sourcearray(NSTEP,NDIM,NGLLX,NGLLY,NGLLZ)

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll


  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
        hgammar(NGLLZ), hpgammar(NGLLZ)
  real(kind=CUSTOM_REAL) :: adj_src(NSTEP,NDIM)

  integer icomp, itime, i, j, k, ios
  double precision :: junk
  character(len=3) :: comp(3)
  character(len=150) :: filename

  call lagrange_any(xi_receiver,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta_receiver,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma_receiver,NGLLZ,zigll,hgammar,hpgammar)

  adj_sourcearray(:,:,:,:,:) = 0.

  comp = (/"BHE", "BHN", "BHZ"/)

  do icomp = 1, NDIM

    filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'
    open(unit = IIN, file = trim(filename), iostat = ios)
    if (ios /= 0) call exit_MPI(myrank, ' file '//trim(filename)//'does not exist')
    do itime = 1, NSTEP
      read(IIN,*) junk, adj_src(itime,icomp)
    enddo
    close(IIN)

  enddo

  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src(:,:)
      enddo
    enddo
  enddo


end subroutine compute_arrays_adjoint_source


