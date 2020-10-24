!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

module signal_processing

  use specfem_par, only: CUSTOM_REAL

  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------------------------------------------------------------------
!>                           butterworth filter
!-----------------------------------------------------------------------------------------------------------------------------------

    subroutine bwfilt (x, y, dt, n, irek, norder, f1, f2)

      ! recursive filtering of data with butterworth filter
      ! x: input array
      ! y: output array
      ! dt: time increment
      ! n: number of data points

      ! irek = 0: forward filtering only
      ! irek = 1: forward and backward filtering

      ! norder: order of butterworth filter
      ! norder = 0: only filtering, no determination of coefficients
      ! norder < 0: no starplots of transfer function and impulse response

      ! f1: low cutoff frequency (Hz)
      ! f1=0: low pass filter

      ! f2: high cutoff frequency (Hz)
      ! f2>0.5/dt: high pass filter

      implicit none
      integer,                                           intent(in)     :: n, norder, irek
      real(kind=CUSTOM_REAL),                            intent(in)     :: dt, f1, f2
      real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout)  :: x, y

      real(kind=CUSTOM_REAL), dimension(10)                             ::  a, b1, b2
      integer                                                           ::  npoles

      if (norder /= 0) then
         npoles=iabs(norder)

         !determination of filter coefficients
         call bpcoeff(f1,f2,npoles, dt, a,b1, b2)

      endif

      if (n /= 0) then
         call rekurs(x,y,n,a,b1,b2,npoles,irek)
      endif

    end subroutine bwfilt

    !---------------------------------------------------------------
    subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)

      ! performs recursive filtering of data in array x of length ndat
      ! filtered output in y
      ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
      ! npoles is the number of poles
      ! iflag=0: forward filtering only
      ! iflag /= 0: forward and backward filtering

      implicit none

      integer,                                           intent(in)    :: ndat, npoles, iflag
      real(kind=CUSTOM_REAL), dimension(10),             intent(in)    :: a, b1, b2
      real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout) :: x, y

      integer                                                          :: n, i
      real(kind=CUSTOM_REAL), dimension(10)                            :: z,z1,z2
      real(kind=CUSTOM_REAL)                                           :: x1,x2


      !forward ---------------------------------------------------------------
      x1 = 0.d0
      x2 = 0.d0

      do i = 1, npoles
         z1(i) = 0.d0
         z2(i) = 0.d0
      enddo

      do n = 1, ndat
         z(1) = a(1)*(x(n)-x2) -b1(1)*z1(1) -b2(1)*z2(1)
         do i = 2, npoles
            z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
         enddo
         x2=x1
         x1=x(n)
         do i = 1, npoles
            z2(i) =z1(i)
            z1(i) =z(i)
         enddo
         y(n) = z(npoles)
      enddo

      if (iflag == 0) then
         return
      endif

      !backward -------------------------------------------------------------
      x1 =0.d0
      x2 =0.d0

      do i = 1, npoles
         z1(i) = 0.d0
         z2(i) = 0.d0
      enddo

      do n = ndat, 1, -1
         z(1) = a(1)*(y(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
         do i =2, npoles
            z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
         enddo
         x2=x1
         x1=y(n)
         do i = 1,npoles
            z2(i)=z1(i)
            z1(i)=z(i)
         enddo
         y(n) = z(npoles)
      enddo

    end subroutine rekurs
    !---------------------------------------------------------------
    subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)

      ! determines filtercoefficients for recursive bandpassfilter

      integer,                                  intent(in)    :: npoles
      real(kind=CUSTOM_REAL),                   intent(in)    :: f1, f2, dt
      real(kind=CUSTOM_REAL),  dimension(10),   intent(inout) :: a,b1,b2

      integer                                                 :: i, npol2, n
      real(kind=CUSTOM_REAL), parameter                       :: pi = 3.141592653589793d0
      real(kind=CUSTOM_REAL)                                  :: d2, w0, w1, w2, ssum, sprod, fact1, fact2, fact3
      complex(kind=CUSTOM_REAL)                               :: s(20), t1, t2, p

      if (npoles > 10) stop 'error: npoles greater than 10'

      d2= 2.d0/dt
      w1=d2*tan(2.d0*pi*f1/d2)
      w2=d2*tan(2.d0*pi*f2/d2)
      w0=0.5*(w2-w1)

      i=1
      npol2=npoles/2+1
      do n =1,npoles
         p = cexp(cmplx(0.d0,dble(2*n-1+npoles)*pi/dble(2*npoles)))
         t1 = p*cmplx(w0,0.d0)
         t2 = sqrt(t1*t1-cmplx(w1*w2,0.d0))
         s(i)=t1+t2
         s(i+1)=t1-t2
         i=i+2
      enddo

      do n=1,npoles
         ssum=2*real(s(n))
         sprod=dble(s(n)*conjg(s(n)))
         fact1=d2*d2-d2*ssum+sprod
         fact2=2.d0*(sprod-d2*d2)
         fact3=d2*d2+d2*ssum+sprod
         a(n)=2.d0*d2*w0/fact1
         b1(n)=fact2/fact1
         b2(n)=fact3/fact1
      enddo

    end subroutine bpcoeff

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------------------------------------------------------------------------------------------------------------
! taper in traces form index i0,..., i3
!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine taper_window_W(t_w,i0,i1,i2,i3,nstep,W)

    implicit none

    integer,                                             intent(in)    :: i0, i1, i2, i3, nstep
    real(kind=CUSTOM_REAL),                              intent(in)    :: W
    real(kind=CUSTOM_REAL), dimension(:), allocatable,   intent(inout) :: t_w

    integer                                                            :: i
    real(kind=CUSTOM_REAL)                                             :: omega, phi, pi

    PI = 3.1415926d0
    t_w(1:nstep)=0.

    ! take off
    omega = pi / (i1 - i0)
    phi = pi / 2 - omega * i1
    do i = i0, i1
       t_w(i) = W*(0.5 + 0.5 *sin(omega * i + phi))
    enddo

    ! flying
    do i = i1+1,i2-1
       t_w(i)=W
    enddo

    ! landing
    omega = pi / (i3 - i2)
    phi = pi/2 - omega * i2
    do i= i2,i3
       t_w(i) = W*(0.5 + 0.5 * sin(omega * i + phi))
    enddo

  end subroutine taper_window_W

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
!> convolution by a wavelet. the output signal have the same length than input
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine convolution_by_wavelet(wavelet, signal, conv_signal, ns, nw)
    implicit none
    integer,                                           intent(in)    :: ns, nw
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(in)    :: wavelet, signal
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout) :: conv_signal
    integer                                                           :: i, k, Ni, Ne

    do i=1, ns
       conv_signal(i)=0._CUSTOM_REAL
       Ni=max(1,i+1-ns)
       Ne=min(i,nw)
       do k=Ni, Ne
          conv_signal(i) = conv_signal(i) + signal(i-k+1) * wavelet(k)
       enddo
    enddo

  end subroutine convolution_by_wavelet
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
!> cross-correlation by a wavelet. the output signal have the same length than input
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine crosscor_by_wavelet(wavelet, signal, conv_signal, ns, nw)
    implicit none
    integer,                                           intent(in)    :: ns, nw
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(in)    :: wavelet, signal
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout) :: conv_signal
    integer                                                           :: i, k, Ni, Ne

    do i=1, ns
       conv_signal(i) = 0._CUSTOM_REAL
       Ni = max(1,i+1-ns)
       Ne = min(i,nw)
       do k = Ni, Ne
          conv_signal(i) = conv_signal(i) + signal(i-k+1) * wavelet(nw - k + 1)
       enddo
    enddo

  end subroutine crosscor_by_wavelet
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
!> compute second time derivative of signal
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine FD2nd(signal, dt, nt)
    implicit none
    real(kind=CUSTOM_REAL),                            intent(in)     :: dt
    integer,                                           intent(in)     :: nt
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout)  :: signal
    real(kind=CUSTOM_REAL), dimension(:), allocatable                 :: wks_signal
    integer                                                           :: i, ier

    allocate(wks_signal(nt),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 296')
    wks_signal(1:nt)=signal(1:nt)

    !! laplacian 1D
    do i=2,nt-1
       signal(i) = ( wks_signal(i+1) - 2._CUSTOM_REAL* wks_signal(i) + wks_signal(i-1) ) / dt / dt
    enddo

    signal(1)=0._CUSTOM_REAL
    signal(nt)=0._CUSTOM_REAL

    deallocate(wks_signal)

  end subroutine FD2nd
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
!> aposidation of signal
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine apodise_sig(signal, nt, lwa)
    implicit none
    integer,                                           intent(in)     :: nt
    real(kind=CUSTOM_REAL),                            intent(in)     :: lwa
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout)  :: signal
    real(kind=CUSTOM_REAL), dimension(:), allocatable                 :: w_tap
    real(kind=CUSTOM_REAL)                                            :: wh
    integer                                                           :: i0, i1, i2, i3, ier

    allocate(w_tap(nt),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 297')
    i0 = 2
    i1 = int(lwa/100.0) * nt + i0
    i3 = nt-1
    i2 = i3  - int(lwa/100.0) * nt
    wh = 1._CUSTOM_REAL
    call taper_window_W(w_tap,i0,i1,i2,i3,nt,wh)
    signal(:) = signal(:)*w_tap(:)

    deallocate(w_tap)

  end subroutine Apodise_sig

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
end module signal_processing
