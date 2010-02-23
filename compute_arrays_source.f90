!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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
        call multiply_arrays_source(sourcearrayd,G11,G12,G13,G21,G22,G23, &
                  G31,G32,G33,hxis,hpxis,hetas,hpetas,hgammas,hpgammas,k,l,m)
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

!================================================================

! we put these multiplications in a separate routine because otherwise
! some compilers try to unroll the six loops above and take forever to compile
  subroutine multiply_arrays_source(sourcearrayd,G11,G12,G13,G21,G22,G23, &
                  G31,G32,G33,hxis,hpxis,hetas,hpetas,hgammas,hpgammas,k,l,m)

  implicit none

  include "constants.h"

! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: G11,G12,G13,G21,G22,G23,G31,G32,G33
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer k,l,m

  integer ir,it,iv

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

  end subroutine multiply_arrays_source

!=============================================================================

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
  real(kind=CUSTOM_REAL),dimension(NSTEP,NDIM,NGLLX,NGLLY,NGLLZ) :: adj_sourcearray

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
        hgammar(NGLLZ), hpgammar(NGLLZ)
        
  real(kind=CUSTOM_REAL) :: adj_src(NSTEP,NDIM)

  integer icomp, itime, i, j, k, ios
  double precision :: junk
  character(len=3),dimension(NDIM) :: comp = (/ "BHN", "BHE", "BHZ" /)
  character(len=256) :: filename

  !adj_sourcearray(:,:,:,:,:) = 0.
  adj_src = 0._CUSTOM_REAL
  
  ! loops over components
  do icomp = 1, NDIM

    filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat = ios)
    if (ios /= 0) cycle ! cycles to next file    
    !if (ios /= 0) call exit_MPI(myrank, ' file '//trim(filename)//'does not exist')
    
    ! reads in adjoint source trace
    do itime = 1, NSTEP
      
      read(IIN,*,iostat=ios) junk, adj_src(itime,icomp)      
      if( ios /= 0 ) &
        call exit_MPI(myrank, &
          'file '//trim(filename)//' has wrong length, please check with your simulation duration')      
    enddo    
    close(IIN)

  enddo

  ! lagrange interpolators for receiver location
  call lagrange_any(xi_receiver,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta_receiver,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma_receiver,NGLLZ,zigll,hgammar,hpgammar)

  ! interpolates adjoint source onto GLL points within this element
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src(:,:)
      enddo
    enddo
  enddo

end subroutine compute_arrays_adjoint_source


! =======================================================================
! compute the integrated derivatives of source parameters (M_jk and X_s)

subroutine compute_adj_source_frechet(displ_s,Mxx,Myy,Mzz,Mxy,Mxz,Myz,eps_s,eps_m_s, &
           hxir,hetar,hgammar,hpxir,hpetar,hpgammar, hprime_xx,hprime_yy,hprime_zz, &
           xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  implicit none

  include 'constants.h'

  ! input
  real(kind=CUSTOM_REAL) :: displ_s(NDIM,NGLLX,NGLLY,NGLLZ)
  double precision :: Mxx, Myy, Mzz, Mxy, Mxz, Myz
  ! output
  real(kind=CUSTOM_REAL) :: eps_s(NDIM,NDIM), eps_m_s(NDIM)

  ! auxilliary
  double precision :: hxir(NGLLX), hetar(NGLLY), hgammar(NGLLZ), &
             hpxir(NGLLX),hpetar(NGLLY),hpgammar(NGLLZ)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

! local variables
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l, tempy1l,tempy2l,tempy3l, &
             tempz1l,tempz2l,tempz3l, hp1, hp2, hp3, &
             xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl, &
             duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl, &
             xix_s,xiy_s,xiz_s,etax_s,etay_s,etaz_s,gammax_s,gammay_s,gammaz_s, &
             hlagrange_xi, hlagrange_eta, hlagrange_gamma, hlagrange

  real(kind=CUSTOM_REAL) :: eps(NDIM,NDIM), eps_array(NDIM,NDIM,NGLLX,NGLLY,NGLLZ), &
             eps_m_array(NGLLX,NGLLY,NGLLZ)

  integer i,j,k,l


! first compute the strain at all the GLL points of the source element
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX

        tempx1l = 0._CUSTOM_REAL
        tempx2l = 0._CUSTOM_REAL
        tempx3l = 0._CUSTOM_REAL

        tempy1l = 0._CUSTOM_REAL
        tempy2l = 0._CUSTOM_REAL
        tempy3l = 0._CUSTOM_REAL

        tempz1l = 0._CUSTOM_REAL
        tempz2l = 0._CUSTOM_REAL
        tempz3l = 0._CUSTOM_REAL

        do l=1,NGLLX
          hp1 = hprime_xx(i,l)
          tempx1l = tempx1l + displ_s(1,l,j,k)*hp1
          tempy1l = tempy1l + displ_s(2,l,j,k)*hp1
          tempz1l = tempz1l + displ_s(3,l,j,k)*hp1

          hp2 = hprime_yy(j,l)
          tempx2l = tempx2l + displ_s(1,i,l,k)*hp2
          tempy2l = tempy2l + displ_s(2,i,l,k)*hp2
          tempz2l = tempz2l + displ_s(3,i,l,k)*hp2

          hp3 = hprime_zz(k,l)
          tempx3l = tempx3l + displ_s(1,i,j,l)*hp3
          tempy3l = tempy3l + displ_s(2,i,j,l)*hp3
          tempz3l = tempz3l + displ_s(3,i,j,l)*hp3
        enddo

! dudx
        xixl = xix(i,j,k)
        xiyl = xiy(i,j,k)
        xizl = xiz(i,j,k)
        etaxl = etax(i,j,k)
        etayl = etay(i,j,k)
        etazl = etaz(i,j,k)
        gammaxl = gammax(i,j,k)
        gammayl = gammay(i,j,k)
        gammazl = gammaz(i,j,k)

        duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
        duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
        duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

        duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
        duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
        duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

        duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
        duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
        duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

! strain eps_jk
        eps(1,1) = duxdxl
        eps(1,2) = (duxdyl + duydxl) / 2
        eps(1,3) = (duxdzl + duzdxl) / 2
        eps(2,2) = duydyl
        eps(2,3) = (duydzl + duzdyl) / 2
        eps(3,3) = duzdzl
        eps(2,1) = eps(1,2)
        eps(3,1) = eps(1,3)
        eps(3,2) = eps(2,3)

        eps_array(:,:,i,j,k) = eps(:,:)

! Mjk eps_jk
        eps_m_array(i,j,k) = Mxx * eps(1,1) + Myy * eps(2,2) + Mzz * eps(3,3) + &
                   2 * (Mxy * eps(1,2) + Mxz * eps(1,3) + Myz * eps(2,3))

      enddo
    enddo
  enddo

  ! interpolate the strain eps_s(:,:) from eps_array(:,:,i,j,k)
  eps_s = 0.
  xix_s = 0.;  xiy_s = 0.;  xiz_s = 0.
  etax_s = 0.; etay_s = 0.; etaz_s = 0.
  gammax_s = 0.; gammay_s = 0.; gammaz_s = 0.

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        hlagrange = hxir(i)*hetar(j)*hgammar(k)

        eps_s(1,1) = eps_s(1,1) + eps_array(1,1,i,j,k)*hlagrange
        eps_s(1,2) = eps_s(1,2) + eps_array(1,2,i,j,k)*hlagrange
        eps_s(1,3) = eps_s(1,3) + eps_array(1,3,i,j,k)*hlagrange
        eps_s(2,2) = eps_s(2,2) + eps_array(2,2,i,j,k)*hlagrange
        eps_s(2,3) = eps_s(2,3) + eps_array(2,3,i,j,k)*hlagrange
        eps_s(3,3) = eps_s(3,3) + eps_array(3,3,i,j,k)*hlagrange

        xix_s = xix_s + xix(i,j,k)*hlagrange
        xiy_s = xiy_s + xiy(i,j,k)*hlagrange
        xiz_s = xiz_s + xiz(i,j,k)*hlagrange
        etax_s = etax_s + etax(i,j,k)*hlagrange
        etay_s = etay_s + etay(i,j,k)*hlagrange
        etaz_s = etaz_s + etaz(i,j,k)*hlagrange
        gammax_s = gammax_s + gammax(i,j,k)*hlagrange
        gammay_s = gammay_s + gammay(i,j,k)*hlagrange
        gammaz_s = gammaz_s + gammaz(i,j,k)*hlagrange

      enddo
    enddo
  enddo

! for completion purpose, not used in specfem3D.f90
  eps_s(2,1) = eps_s(1,2)
  eps_s(3,1) = eps_s(1,3)
  eps_s(3,2) = eps_s(2,3)

! compute the gradient of M_jk * eps_jk, and then interpolate it

  eps_m_s = 0.
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        hlagrange_xi = hpxir(i)*hetar(j)*hgammar(k)
        hlagrange_eta = hxir(i)*hpetar(j)*hgammar(k)
        hlagrange_gamma = hxir(i)*hetar(j)*hpgammar(k)

        eps_m_s(1) = eps_m_s(1) +  eps_m_array(i,j,k) * (hlagrange_xi * xix_s &
                   + hlagrange_eta * etax_s + hlagrange_gamma * gammax_s)
        eps_m_s(2) = eps_m_s(2) +  eps_m_array(i,j,k) * (hlagrange_xi * xiy_s &
                   + hlagrange_eta * etay_s + hlagrange_gamma * gammay_s)
        eps_m_s(3) = eps_m_s(3) +  eps_m_array(i,j,k) * (hlagrange_xi * xiz_s &
                   + hlagrange_eta * etaz_s + hlagrange_gamma * gammaz_s)

      enddo
    enddo
  enddo

end subroutine compute_adj_source_frechet

! =======================================================================

! compute array for acoustic source
  subroutine compute_arrays_source_acoustic(xi_source,eta_source,gamma_source,&
                        sourcearray,xigll,yigll,zigll,factor_source)

  implicit none

  include "constants.h"

  double precision :: xi_source,eta_source,gamma_source
  real(kind=CUSTOM_REAL) :: factor_source
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

! local parameters
! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas
  integer :: i,j,k
  
! initializes  
  sourcearray(:,:,:,:) = 0._CUSTOM_REAL
  sourcearrayd(:,:,:,:) = 0.d0

! computes Lagrange polynomials at the source location
  call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  call lagrange_any(eta_source,NGLLY,yigll,hetas,hpetas)
  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)

! calculates source array for interpolated location
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        ! identical source array components in x,y,z-direction
        sourcearrayd(:,i,j,k) = hxis(i)*hetas(j)*hgammas(k)*dble(factor_source)        
      enddo
    enddo
  enddo

! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then
    sourcearray(:,:,:,:) = sngl(sourcearrayd(:,:,:,:))
  else
    sourcearray(:,:,:,:) = sourcearrayd(:,:,:,:)
  endif

  end subroutine compute_arrays_source_acoustic


