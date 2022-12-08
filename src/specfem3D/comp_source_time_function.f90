!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

  double precision function comp_source_time_function(t,hdur)

  use specfem_par, only: DT,t0,PI

  implicit none

  double precision, intent(in) :: t,hdur

  double precision, external :: netlib_specfun_erf

  ! taper - experimental test
  ! idea: to see if tapering the onset of the source time function would diminish numerical noise.
  !       so far not really, spurious oscillations still occur due to discretization.
  !       thus, artefacts seem not too much affected by non-zero onset values of source time function.
  logical, parameter :: USE_TAPERED_BEGINNING = .false.
  double precision, parameter :: length_window = 200
  double precision :: taper
  integer :: l

  ! quasi Heaviside, small Gaussian moment-rate tensor with hdur
  comp_source_time_function = 0.5d0*(1.0d0 + netlib_specfun_erf(t/hdur))

  ! taper to suppress noise?
  if (USE_TAPERED_BEGINNING) then
    if (t < length_window * DT - t0) then
      ! check if taper length is reasonable compared to hdur
      if ( nint(hdur/DT) > length_window) then
        ! index from length_window to 0
        l = nint( ((length_window * DT - t0) - t)/DT )
        ! from 0 to length_window-1
        l = length_window - l
        ! cosine taper, otherwise using a constant (1.0) instead
        taper = (1.d0 - cos(PI*(l+1)/(length_window)))/2.d0
        ! linear taper
        !taper = dble(l+1)/dble(length_window)
        ! tapers stf
        comp_source_time_function = comp_source_time_function * taper
        !debug
        !print *,'debug: taper ',taper,l,t,length_window*DT-t0,hdur,nint(hdur/DT)
      endif
    endif
  endif

  end function comp_source_time_function


!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_gauss(t,hdur)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,hdur
  double precision :: hdur_decay,a

  ! note: hdur given is hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
  !           and SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.68
  hdur_decay = hdur

  ! this here uses a stronger Gaussian decay rate (empirical value) to avoid non-zero onset times;
  ! however, it should mimik a triangle source time function...
  !hdur_decay = hdur  / SOURCE_DECAY_STRONG

  ! note: a nonzero time to start the simulation with would lead to more high-frequency noise
  !          due to the (spatial) discretization of the point source on the mesh

  ! Gaussian wavelet
  a = 1.d0 / (hdur_decay**2)
  comp_source_time_function_gauss = exp(-a*t**2) / (sqrt(PI)*hdur_decay)

  end function comp_source_time_function_gauss

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_gauss_2(t,hdur)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,hdur

  ! local parameters
  double precision :: a

  ! source time function as defined in:
  ! M.A. Meschede, C.L. Myhrvold and J. Tromp, 2011.
  ! Antipodal focusing of seismic waves due to large meteorite impacts on Earth,
  ! GJI, 187, p. 529-537
  !
  ! equation (2):
  ! S(t) = sqrt(pi/tau**2) * exp(-pi**2 * t**2 / tau**2)

  ! factor
  a = sqrt(PI / hdur**2)

  comp_source_time_function_gauss_2 = a * exp(-PI**2 * t**2 / hdur**2)

  end function comp_source_time_function_gauss_2

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_dgau(t,hdur)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,hdur
  double precision :: hdur_decay,a

  ! note: hdur given is hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
  !           and SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.68
  hdur_decay = hdur

  ! this here uses a stronger Gaussian decay rate (empirical value) to avoid non-zero onset times;
  ! however, it should mimik a triangle source time function...
  !hdur_decay = hdur  / SOURCE_DECAY_STRONG

  ! note: a nonzero time to start the simulation with would lead to more high-frequency noise
  !          due to the (spatial) discretization of the point source on the mesh

  ! first derivative of a Gaussian wavelet
  a = 1.d0 / (hdur_decay**2)
  comp_source_time_function_dgau = - 2.d0 * a * t * exp(-a*t**2) / (sqrt(PI)*hdur_decay)

  end function comp_source_time_function_dgau

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_d2gau(t,hdur)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,hdur
  double precision :: hdur_decay,a

  ! note: hdur given is hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
  !           and SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.68
  hdur_decay = hdur

  ! this here uses a stronger Gaussian decay rate (empirical value) to avoid non-zero onset times;
  ! however, it should mimik a triangle source time function...
  !hdur_decay = hdur  / SOURCE_DECAY_STRONG

  ! note: a nonzero time to start the simulation with would lead to more high-frequency noise
  !          due to the (spatial) discretization of the point source on the mesh

  ! second derivative of a Gaussian wavelet
  a = 1.d0 / (hdur_decay**2)
  comp_source_time_function_d2gau = 2.d0 * a * (-1.d0 + 2.d0*a*t**2) * exp(-a*t**2) / (sqrt(PI)*hdur_decay)

  end function comp_source_time_function_d2gau

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_rickr(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0

  ! local variables
  double precision :: a

  ! Ricker wavelet
  a = PI**2 * f0**2
  comp_source_time_function_rickr = (1.d0 - 2.d0*a*t*t) * exp( -a*t*t )

  !!! another source time function they have called 'ricker' in some old papers,
  !!! e.g., 'Finite-Frequency Kernels Based on Adjoint Methods' by Liu & Tromp, BSSA (2006)
  !!! in order to benchmark those simulations, the following formula is needed.
  ! comp_source_time_function_rickr = -2.d0*PI*PI*f0*f0*f0*t * exp(-PI*PI*f0*f0*t*t)

  end function comp_source_time_function_rickr

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_drck(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0

  ! local variables
  double precision :: a

  ! first derivative of a Ricker wavelet
  a = PI**2 * f0**2
  comp_source_time_function_drck = 2.d0*a*t * (-3.d0 + 2.d0*a*t*t) * exp( -a*t*t )

  end function comp_source_time_function_drck

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_d2rck(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0

  ! local variables
  double precision :: a

  ! second derivative of a Ricker wavelet
  a = PI**2 * f0**2
  comp_source_time_function_d2rck = -2.d0*a * (3.d0 - 12.d0*a*t*t + 4.d0*a**2*t*t*t*t) * exp( -a*t*t )

  end function comp_source_time_function_d2rck

!
!-------------------------------------------------------------------------------------------------
!


  double precision function comp_source_time_function_mono(t,f0)

  use constants, only: PI,TAPER_MONOCHROMATIC_SOURCE

  implicit none

  double precision,intent(in) :: t,f0
  double precision :: tt
  integer :: taper

  tt = 2 * PI * f0 * t
  taper = ceiling(TAPER_MONOCHROMATIC_SOURCE * f0)

  if (t < taper / f0) then
    comp_source_time_function_mono = sin(tt) * (0.5 - 0.5 * cos(tt / taper / 2.0))
  else
    comp_source_time_function_mono = sin(tt)
  endif

  ! monochromatic source time function

  end function comp_source_time_function_mono

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_ext(it_index,isource)

  use specfem_par, only: myrank, NSTEP, user_source_time_function

  implicit none

  integer,intent(in) :: it_index,isource

  ! safety check
  if (.not. allocated (user_source_time_function)) then
    stop 'Error invalid user STF array, not allocated yet'
  endif

  ! checks index bounds
  if (it_index < 1 .or. it_index > NSTEP) then
    print *,'Error: external source time function index ',it_index,'should be between 1 and ',NSTEP
    call exit_MPI(myrank,'Invalid external source time function index in comp_source_time_function_ext() routine')
  endif

  ! gets stored STF
  comp_source_time_function_ext = user_source_time_function(it_index,isource)

  end function comp_source_time_function_ext
