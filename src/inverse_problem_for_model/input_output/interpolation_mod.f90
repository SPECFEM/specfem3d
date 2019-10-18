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

module interpolation_mod

  use constants, only: CUSTOM_REAL

  implicit none

  !**************************************************
  !* CONTENT OF MY MODULE PRECISION_MOD
  !*** Single precision
  integer, parameter :: si = selected_int_kind(8)
  integer, parameter :: sp = selected_real_kind(4)

  !*** Double precision
  integer, parameter :: di = selected_int_kind(16)
  integer, parameter :: dp = selected_real_kind(8)

  !*** Custom precision
!!!!!!!! DK DK  integer, parameter :: cp = selected_real_kind(4)
  integer, parameter :: cp = CUSTOM_REAL

  !*** Special precision
  integer, parameter :: hp = selected_real_kind(8)
  !***************************************************

  !***************************************************
  !* CONTENT OF MY MODULE CONSTANTS_MOD
  !*** Constants for inititalization
  real(kind=cp), parameter :: mytiny = 1.e-9_cp
  real(kind=cp), parameter :: myhuge = 1.e30_cp

  real(kind=cp),    parameter :: myverysmall = 1.e-25_cp
  real(kind=dp),    parameter :: mysmall     = 0.000001_dp
  integer(kind=si), parameter :: myhugeint   = 100000000

  real(kind=cp), parameter :: zero = 0._cp
  real(kind=cp), parameter :: one  = 1._cp
  real(kind=cp), parameter :: half = 0.5_cp

  real(kind=cp), parameter :: onethird  = 1._cp/3._cp
  real(kind=cp), parameter :: twothird  = 2._cp/3._cp
  real(kind=cp), parameter :: fourthird = 4._cp/3._cp

  real(kind=cp), parameter :: onefourth = 1._cp/4._cp
  real(kind=cp), parameter :: onesixth  = 1._cp/6._cp

  !*** Constants
  real(kind=cp), parameter :: pi  = 3.141592653589793_cp
  real(kind=cp), parameter :: twopi = 2._cp * pi

  !*** Constants for conversions
  real(kind=cp), parameter :: deg2rad = 0.017453292519943
  real(kind=cp), parameter :: rad2deg = 57.295779513082323

  !*** Check stability in modeling
  real(kind=cp), parameter :: stability_criterion = 1.e28_cp
  !**********************************************************


contains


!================================================================================
! Filtering routines
  subroutine bwfilt (x, y, dt, n, irek, norder, f1, f2)

    ! recursive filtering of data with butterworth filter
    ! x: input array
    ! y: output array
    ! dt: time increment
    ! n: number of data points

    ! irek=0: forward filtering only
    ! irek=1: forward and backward filtering

    ! norder: order of butterworth filter
    ! norder=0: only filtering, no determination of coefficients
    ! norder < 0: no starplots of transfer function and impulse response

    ! f1: low cutoff frequency (Hz)
    ! f1=0: low pass filter

    ! f2: high cutoff frequency (Hz)
    ! f2>0.5/dt: high pass filter

    implicit none

    integer(kind=si) :: iunit, npoles,n,lx
    integer(kind=si) :: irek,norder
    real(kind=cp), dimension(n) ::x,y
    real(kind=cp), dimension (10) ::  a, b1, b2
    real(kind=cp) :: dt,f1,f2

    !real(kind(0d0)) :: x(n),y(n)

    iunit = 3

    if (norder /= 0) then
       npoles=abs(norder)
       !determination of filter coefficients
       call bpcoeff(f1,f2,npoles, dt, a,b1, b2)
       if (norder >= 0) then
          !plot of transfer function and impuulse response
          lx = 100
          !filtering
       endif
    endif


    if (n /= 0) then
       call rekurs(x,y,n,a,b1,b2,npoles,irek)
    endif
    return
  end subroutine bwfilt

  subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
    ! performs recursive filtering of data in array x of length ndat
    ! filtered output in y
    ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
    ! npoles is the number of poles
    ! iflag=0: forward filtering only
    ! iflag /= 0: forward and backward filtering

    implicit none

    real(kind=cp), dimension(10) :: z,z1,z2 ,a,b1,b2
    real(kind=cp)  ::  x1,x2
    integer(kind=si) :: ndat, npoles,n,i
    integer(kind=si) :: iflag
    real(kind=cp), dimension(ndat) :: x, y

    !forward

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

    !backward

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
    return
  end subroutine rekurs

  subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
    !determines filtercoefficients for recursive bandpassfilter

    real(kind=cp),dimension(10) :: a,b1,b2
    complex(kind=4) :: s(20), t1,t2,p
    real(kind=cp), parameter :: pi = 3.141592653589793
    real(kind=cp) :: f1,f2,dt,d2,w0,w1,w2,ssum, sprod,fact1,fact2,fact3
    integer(kind=si) :: i,npol2,n,npoles


    if (npoles > 10) then
       stop ' npoles greater than 10: STOP '
    endif

    d2= 2/dt
    w1=d2*tan(2.*pi*f1/d2)
    w2=d2*tan(2.*pi*f2/d2)
    w0=0.5*(w2-w1)

    i=1
    npol2=npoles/2+1
    do n =1,npoles
       p = cexp(cmplx(0.,real(2*n-1+npoles)*pi/real(2*npoles)))
       t1 = p*cmplx(w0,0.)
       t2 = sqrt(t1*t1-cmplx(w1*w2,0.))
       s(i)=t1+t2
       s(i+1)=t1-t2
       i=i+2
    enddo

    do n=1,npoles
       ssum=2*real(s(n))
       sprod=real(s(n)*conjg(s(n)))
       fact1=d2*d2-d2*ssum+sprod
       fact2=2.*(sprod-d2*d2)
       fact3=d2*d2+d2*ssum+sprod
       a(n)=2.*d2*w0/fact1
       b1(n)=fact2/fact1
       b2(n)=fact3/fact1
    enddo
    return
  end subroutine bpcoeff
!--------------------------------------------------------------------------------

!================================================================================
! Iterative time domain deconvolution
  subroutine time_deconv(dobs,dcal,dt,nt,nit,src_sum)

    integer(kind=si), intent(in) :: nt, nit
    integer(kind=si)             :: i, ii
    integer(kind=si), parameter  :: part=0

    real(kind=cp), intent(in) :: dt
    real(kind=cp)             :: amp_corr

    real(kind=cp), dimension(nt), intent(in)    :: dobs, dcal
    real(kind=cp), dimension(:), allocatable, intent(inout) :: src_sum

    real(kind=cp), dimension(:), allocatable :: dobs2, autocorr, crosscorr, new_obs, src_one

    integer :: ier

    !* 0. Allocate
    if (.not. allocated(autocorr)) then
      allocate(autocorr(nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 377')
    endif
    if (.not. allocated(crosscorr)) then
      allocate(crosscorr(nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 378')
    endif
    if (.not. allocated(src_sum)) then
      allocate(src_sum(nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 379')
    endif
    if (.not. allocated(src_one)) then
      allocate(src_one(nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 380')
    endif
    if (.not. allocated(dobs2)) then
      allocate(dobs2(nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 381')
    endif
    if (.not. allocated(new_obs)) then
      allocate(new_obs(nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 382')
    endif

    src_sum = 0._cp
    src_one = 0._cp

    !* 1. Estimate auto-corr of computed data
    call mycorrelation(dcal,dcal,nt,nt,autocorr,part)
    autocorr = autocorr * dt
    amp_corr = 1._cp / maxval(abs(autocorr))

    !* 2. Time iteration to estimate sources
    dobs2 = dobs
    do i = 1, nit

       !* Compute correlation
       call mycorrelation(dobs2,dcal,nt,nt,crosscorr,part)
       crosscorr = crosscorr * dt

       !* Find maximum of correlation
       ii = maxloc(abs(crosscorr),dim=1)

       !* Put local contibution to src_one and src_sum
       src_one     = 0._cp
       src_one(ii) = crosscorr(ii) * amp_corr
       src_sum     = src_sum + src_one

       !* Convolve
       call myconvolution(src_one,dcal,nt,nt,new_obs,part)
       new_obs = new_obs * dt
       dobs2 = dobs2 - new_obs

    enddo

  end subroutine time_deconv

!--------------------------------------------------------------------------------

!================================================================================
! Convolution routine
   subroutine myconvolution(sig1,sig2,n1,n2,conv,part)

     integer(kind=si), intent(in) :: n1, n2, part

     real(kind=cp), dimension(n1), intent(in) :: sig1
     real(kind=cp), dimension(n2), intent(in) :: sig2

     real(kind=cp), dimension(:), allocatable, intent(inout) ::conv

     real(kind=cp), dimension(n1+n2-1) :: convtmp

     integer(kind=si) :: i1, i2, ind

     integer :: ier

     !*** Put to zero
     convtmp = zero

     !*** Convolve
     do i1=1,n1
        do i2=1,n2
           convtmp(i1+i2-1) = convtmp(i1+i2-1) + sig1(i1) * sig2(i2)
        enddo
     enddo

     if (part == 0) then !(middle regular)
        !*** Take good parts (this is wrong...)
        if (modulo(n2,2) == 0) then
           ind = n2/2+1
        else
           ind = ceiling(real(n2/2,kind=cp))
        endif
        if (.not. allocated(conv)) then
          allocate(conv(n2),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 383')
        endif
        conv(1:n2) = convtmp(ind:ind+n2-1)

     else if (part == 1) then ! full convolution

        if (.not. allocated(conv)) then
          allocate(conv(n1+n2-1),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 384')
        endif
        conv(:) = convtmp(:)

     else if (part == 2) then !(middle irregular)
        if (.not. allocated(conv)) then
          allocate(conv(n2),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 385')
        endif
        conv(1:n2) = convtmp(n2:n2+n2-1)
     endif

   end subroutine myconvolution
!--------------------------------------------------------------------------------

!================================================================================
! Correlation routine
   subroutine mycorrelation(sig1,sig2,n1,n2,corr,part)

     integer(kind=si), intent(in) :: n1, n2, part

     real(kind=cp), dimension(n1), intent(in) :: sig1
     real(kind=cp), dimension(n2), intent(in) :: sig2

     real(kind=cp), dimension(:), allocatable, intent(inout) :: corr

     real(kind=cp), dimension(n2) :: flipsig2
     integer(kind=si) :: i

     integer :: ier

     !*** Choose size of corr
     if (part == 0) then ! (middle)
        if (.not. allocated(corr)) then
          allocate(corr(n2),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 386')
        endif
     else if (part == 1) then ! (full)
        if (.not. allocated(corr)) then
          allocate(corr(n1+n2-1),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 387')
        endif
     endif

     !*** Flip second signal
     do i=1,n2
        flipsig2(i) = sig2(n2-i+1)
     enddo

     !*** Use regular convolution
     call myconvolution(sig1,flipsig2,n1,n2,corr,part)

   end subroutine mycorrelation
!--------------------------------------------------------------------------------

!================================================================================
! Find lag
   subroutine determine_lag(sig1,sig2,n1,n2,lag)

     integer(kind=si) :: part=0

     integer(kind=si), intent(in) :: n1, n2

     real(kind=cp), dimension(n1), intent(in) :: sig1
     real(kind=cp), dimension(n2), intent(in) :: sig2

     real(kind=cp), dimension(:), allocatable :: corr

     integer(kind=si), intent(out) :: lag

     integer(kind=si) :: it, ind
     real(kind=cp)    :: maxcorr

     integer :: ier

     !*** Take good parts, define middle
     if (.not. allocated(corr)) then
       allocate(corr(n1),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 388')
     endif

     if (modulo(n2,2) == 0) then
        ind = n2/2
     else
        ind = (n2-1)/2 +1
     endif

     !*** Cross-correlation
     call mycorrelation(sig1,sig2,n1,n2,corr,part)

     !*** Find delay
     maxcorr = maxval(abs(corr))
     do it=1,n1
        if (abs(corr(it)) == maxcorr) then
           lag = it-ind
           exit
!        else
!           stop 'there is a problem....'
        endif
     enddo

   end subroutine determine_lag
!--------------------------------------------------------------------------------

!================================================================================
! Compute stalta ratio and give first pick
   subroutine substalta(sig,n,nsta,nlta,crit,stalta,tpick)

     integer(kind=si), intent(in) :: n, nsta, nlta
     real(kind=cp), intent(in) :: crit
     real(kind=cp), dimension(n), intent(in) :: sig

     integer(kind=si), intent(out) :: tpick

     real(kind=cp), dimension(n) :: sta, lta
     real(kind=cp), dimension(n), intent(out) :: stalta
     real(kind=cp), dimension(n+2*nsta) :: tmpsta
     real(kind=cp), dimension(n+2*nlta) :: tmplta

     integer(kind=si) :: i

     !*** Compute the short time average (STA)
     tmpsta(1:nsta) = sig(1)
     tmpsta(nsta+1:nsta+n) = sig(:)
     tmpsta(nsta+n+1:n+2*nsta) = sig(n)
     sta = zero
     do i=1+nsta,n+nsta
        sta(i-nsta) = sum(tmpsta(i-nsta:i+nsta)**2)
     enddo
     sta = 0.5 * sta / nsta

     !*** Compute the long time average (LTA)
     tmplta(1:nlta) = sig(1)
     tmplta(nlta+1:nlta+n) = sig(:)
     tmplta(nlta+n+1:n+2*nlta) = sig(n)
     lta = zero
     do i=1+nlta,n+nlta
        lta(i-nlta) = sum(tmplta(i-nlta:i+nlta)**2)
     enddo
     lta = 0.5 * lta / nlta

     !*** Compute ratio and gives first pick
     stalta = sta / lta
     do i=1,n
        if (stalta(i) >= crit) then
           tpick = i
           exit
        else
           tpick = n
        endif
     enddo

   end subroutine substalta
!--------------------------------------------------------------------------------

!================================================================================
! Determine first arrival
   subroutine pick_first_arrival(sig1,n1,tpick,dt)

     integer(kind=si), intent(in) :: n1
     real(kind=cp)   , intent(in) :: dt
     real(kind=cp), dimension(n1), intent(in) :: sig1

     integer(kind=si), intent(out) :: tpick

     real(kind=cp), dimension(n1) :: stalta

     integer(kind=si) :: nsta, nlta
     real(kind=cp)    :: crit

     nsta = int(ceiling(5. / dt))   !* 5s pour sta
     nlta = int(ceiling(60. / dt))  !* 60s pour lta

     crit = 2.1                     !* a bit more than twice...

     call substalta(sig1,n1,nsta,nlta,crit,stalta,tpick)

   end subroutine pick_first_arrival
!--------------------------------------------------------------------------------

!================================================================================

! Trilinear interpolation

  subroutine trilin_interp(x,y,z,valx,valy,valz,nx,ny,nz,valin,valout)!,indx,indy,indz)

    !*** Size of grid
    integer(kind=si), intent(in) :: nx, ny, nz

    !*** Lookup 1D tables with Cartesian coordinates
    real(kind=cp), dimension(nx), intent(in) :: valx
    real(kind=cp), dimension(ny), intent(in) :: valy
    real(kind=cp), dimension(nz), intent(in) :: valz

    !*** Input grid of data
    real(kind=cp), dimension(nx,ny,nz), intent(in) :: valin

    !*** Coordinate for interpolated value
    real(kind=cp), intent(in) :: x, y, z

    !*** Temporary variables
    real(kind=cp) :: xix1, xix2, yiy1, yiy2, ziz1, ziz2
    real(kind=cp) :: x2x1, x1x2, y2y1, y1y2, z2z1, z1z2
    real(kind=cp) :: facx1y1z1, facx1y1z2, facx1y2z1, facx1y2z2, facx2y1z1, facx2y1z2, facx2y2z1, facx2y2z2

    !*** Output value to be interpolated
    real(kind=cp), intent(out) :: valout
    integer(kind=si) :: indx, indy, indz  !!! Previous guest
    integer(kind=si) :: kx, ky, kz, m

    !*** Lookup table
!    call locate_hunt(valx,nx,x,indx)
!    call locate_hunt(valy,ny,y,indy)
!    call locate_hunt(valz,nz,z,indz)
    call locate_bissection(valx,nx,x,indx)
    call locate_bissection(valy,ny,y,indy)
    call locate_bissection(valz,nz,z,indz)


    m=2
    kx = min(max(indx-(m-1)/2,1),nx+1-m)
    ky = min(max(indy-(m-1)/2,1),ny+1-m)
    kz = min(max(indz-(m-1)/2,1),nz+1-m)


    !*** x_i - x1
    xix1 = x - valx(kx)
    yiy1 = y - valy(ky)
    ziz1 = z - valz(kz)

    !*** x_i - x2
    xix2 = x - valx(kx+1)
    yiy2 = y - valy(ky+1)
    ziz2 = z - valz(kz+1)

    !*** x1 - x2
    x1x2 = 1./ (valx(kx) - valx(kx+1))
    y1y2 = 1./ (valy(ky) - valy(ky+1))
    z1z2 = 1./ (valz(kz) - valz(kz+1))

    !*** x2 - x1
    x2x1 = 1./(valx(kx+1) - valx(kx))
    y2y1 = 1./(valy(ky+1) - valy(ky))
    z2z1 = 1./(valz(kz+1) - valz(kz))

    !*** Factors
    facx1y1z1 = xix2*yiy2*ziz2 * x1x2*y1y2*z1z2 * valin(kx,ky,kz)
    facx1y1z2 = xix2*yiy2*ziz1 * x1x2*y1y2*z2z1 * valin(kx,ky,kz+1)
    facx1y2z1 = xix2*yiy1*ziz2 * x1x2*y2y1*z1z2 * valin(kx,ky+1,kz)
    facx1y2z2 = xix2*yiy1*ziz1 * x1x2*y2y1*z2z1 * valin(kx,ky+1,kz+1)
    facx2y1z1 = xix1*yiy2*ziz2 * x2x1*y1y2*z1z2 * valin(kx+1,ky,kz)
    facx2y1z2 = xix1*yiy2*ziz1 * x2x1*y1y2*z2z1 * valin(kx+1,ky,kz+1)
    facx2y2z1 = xix1*yiy1*ziz2 * x2x1*y2y1*z1z2 * valin(kx+1,ky+1,kz)
    facx2y2z2 = xix1*yiy1*ziz1 * x2x1*y2y1*z2z1 * valin(kx+1,ky+1,kz+1)

    !*** Final value
    valout = facx1y1z1 + facx1y1z2 + facx1y2z1 + facx1y2z2 + facx2y1z1 + facx2y1z2 + facx2y2z1 + facx2y2z2

  end subroutine trilin_interp
!--------------------------------------------------------------------------------



!================================================================================
! Locate with bisection
  subroutine locate_bissection(valx,n,x,ind)

    integer(kind=si), intent(in) :: n
    real(kind=cp),    intent(in) :: x

    real(kind=cp), dimension(n), intent(in) :: valx

    integer(kind=si), intent(out) :: ind

    integer(kind=si) :: jl, ju, jm

    jl = 0
    ju = n+1

    do while ( (ju-jl) > 1 )
       jm = (ju + jl) / 2
       if ( x >= valx(jm) ) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if ( x == valx(1) ) then
       ind = 1
    else if ( x == valx(n) ) then
       ind = n-1
    else
       ind = jl
    endif

  end subroutine locate_bissection
!--------------------------------------------------------------------------------


!================================================================================
! Locate with hunt
  subroutine locate_hunt(valx,n,x,jl)

    integer(kind=si), intent(in) :: n
    real(kind=cp),    intent(in) :: x

    real(kind=cp), dimension(n), intent(in) :: valx

    integer(kind=si), intent(inout) :: jl    ! inout because initial guess

    integer(kind=si) :: ju, jm,  inc

    inc = 1                      ! Set the hunting increment

    if (jl <= 0 .or. jl > n) then ! input guess not useful goto bissection
       jl = 0
       ju = n+1
    else
       if (x >= valx(jl)) then         ! Hunt up
          do while ( x >= valx(jl) )   ! Not done with hunting
             ju  = jl + inc
             if ( ju > n ) then        ! Done hunting... end of table
                ju = n+1
                exit
             else if (x < valx(ju)) then    ! Found bracket
                exit
             else                      ! Not done => double increment
                jl = ju
                inc = inc + inc
             endif
          enddo
       else                            ! Hunt down
          ju = jl
          do while ( x < valx(jl) )    ! not done hunting
             jl = jl - inc
             if (jl <= 0) then         ! Done hunting... end of table
                jl= 0
                exit
             else if (x >= valx(ju)) then   ! Found bracket
                exit
             else
                ju = jl
                inc = inc + inc
             endif
          enddo
       endif
    endif

    !** hunt is done
    ! Dichotomy
    do while ( (ju-jl) > 1 )
       jm = (ju + jl) / 2
       if ( x >= valx(jm) ) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if ( x == valx(1) ) then
       jl = 1
    else if ( x == valx(n) ) then
       jl = n-1
    endif

  end subroutine locate_hunt
!--------------------------------------------------------------------------------


!================================================================================
! Taper 3D
  subroutine taper_3D(ndom,taper,isr,sizetapx,sizetapy,sizetapz)

     integer(kind=si), dimension(3), intent(in) :: ndom

     integer(kind=si) :: i, j, k, isr

     real(kind=cp), dimension(:,:,:), allocatable, intent(inout) :: taper

     real(kind=cp), dimension(:), allocatable :: tapx, tapy, tapz
     real(kind=cp) :: alpha
     real(kind=cp),intent(in) :: sizetapx, sizetapy, sizetapz

     integer :: ier

     if (.not. allocated(tapx)) then
       allocate(tapx(ndom(1)),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 389')
     endif
     if (.not. allocated(tapy)) then
       allocate(tapy(ndom(2)),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 390')
     endif
     if (.not. allocated(tapz)) then
       allocate(tapz(ndom(3)),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 391')
     endif

     alpha = sizetapx*2./ndom(1)
     tapx = tuckeywin(ndom(1),alpha)

     alpha = sizetapy*2./ndom(2)
     tapy = tuckeywin(ndom(2),alpha)

     alpha = sizetapz*2./ndom(3)
     tapz = tuckeywin(ndom(3),alpha)

     !*** Remove taper at the top
     if (isr == 2) then
        tapz(ndom(3)/2:ndom(3)) = 1.
        tapz(ndom(3))   = 0.1
        tapz(ndom(3)-1) = 0.2
        tapz(ndom(3)-2) = 0.3
        tapz(ndom(3)-3) = 0.4
        tapz(ndom(3)-4) = 0.5
        tapz(ndom(3)-5) = 0.6
        tapz(ndom(3)-6) = 0.7
        tapz(ndom(3)-7) = 0.8
        tapz(ndom(3)-8) = 0.9

     else if (isr == 1) then
                tapz(ndom(3)/2:ndom(3))=1.
     endif

     do k=1,ndom(3)
        do j=1,ndom(2)
           do i=1,ndom(1)
              taper(i,j,k) = tapx(i) * tapy(j) * tapz(k)
           enddo
        enddo
     enddo

     if (allocated(tapx)) deallocate(tapx)
     if (allocated(tapy)) deallocate(tapy)
     if (allocated(tapz)) deallocate(tapz)

  end subroutine taper_3D
!--------------------------------------------------------------------------------


!================================================================================
! Tukey tapering windows
!--------------------------------------------------
! N     : number of samples
! alpha : percentage of signal to taper
! tuk   : tapered window
  function tuckeywin(N,alpha) result(tuk)

    integer(kind=si), intent(in)   :: N
    real(kind=cp), intent(in)      :: alpha

    integer(kind=si) :: i
    real(kind=cp), parameter :: pipi=3.141592653589793
    real(kind=cp), dimension(N) :: tuk

    !*** Central part
    tuk(:) = 1.

    !*** Left part
    do i=0,int(0.5*alpha*(N-1))
       if (i+1 > 0 .and. i+1 <= N) then
          tuk(i+1) = 0.5*(1+cos(pipi*(2.*i/(alpha*(N-1.))-1.)))
       endif
    enddo

    !*** Right part
    do i=int((N-1)*(1-alpha/2.)),N-1
       if (i+1 > 0 .and. i+1 <= N) then
          tuk(i+1) = 0.5*(1+cos(pipi*(2.*i/(alpha*(N-1.))-(2./alpha)+1.)))
       endif
    enddo

  end function tuckeywin
!--------------------------------------------------------------------------------



!================================================================================
! Cardinal sine function
!--------------------------------------------------
! x = argument of sinc function
  real(kind=dp) function mysinc(x)

    real(kind=dp) :: x
    real(kind=dp), parameter :: pipi=3.141592653589793

    if (abs(x) >= 1e-13) then
       mysinc = sin(pipi*x)/(pipi*x)
    else
       mysinc = 1._dp
    endif

  end function mysinc
!--------------------------------------------------------------------------------

!!$!================================================================================
!!$! Make lookup tables
!!$  subroutine create_lookup_tables(nx,ny,nz,xo,yo,zo,dx,dy,dz, &
!!$                                  xcrd, ycrd, zcrd)
!!$
!!$    integer(kind=si), intent(in) :: nx, ny, nz
!!$    real(kind=cp),    intent(in) :: xo, yo, zo
!!$    real(kind=cp),    intent(in) :: dx, dy, dz
!!$
!!$    integer(kind=si) :: xmid, ymid, zmid, i
!!$
!!$    real(kind=cp), dimension(nx), intent(out) :: xcrd
!!$    real(kind=cp), dimension(ny), intent(out) :: ycrd
!!$    real(kind=cp), dimension(nz), intent(out) :: zcrd
!!$
!!$    !*** Creates coordinates
!!$    !* For x
!!$    if (modulo(nx,2)==0) then
!!$       xmid = nx / 2
!!$    else
!!$       xmid = (nx / 2) + 1
!!$    endif
!!$
!!$    do i=1,nx
!!$       xcrd(i) = -dble(xmid - i) * dx
!!$    enddo
!!$
!!$    !* For y
!!$    if (modulo(ny,2)==0) then
!!$       ymid = ny / 2
!!$    else
!!$       ymid = (ny / 2) + 1
!!$    endif
!!$
!!$    do i=1,ny
!!$       ycrd(i) = -dble(ymid - i) * dy
!!$    enddo
!!$
!!$    !* For z
!!$    zmid = 1       ! Indice of zorigine = 0
!!$    do i=1,nz
!!$       zcrd(i) = -model_size(3) + dble(i-zmid) * inv_dz
!!$    enddo
!!$
!!$  end subroutine create_inversion_grid

end module interpolation_mod
