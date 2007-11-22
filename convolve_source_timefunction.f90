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

  program convolve_source_time_function

!
! convolve seismograms computed for a Heaviside with given source time function
!

! we mimic a triangle of half duration equal to half_duration_triangle
! using a Gaussian having a very close shape, as explained in Figure 4.2
! of the manual

  implicit none

  include "constants.h"

  integer :: i,j,N_j,number_remove,nlines

  double precision :: alpha,dt,tau_j,source,exponent,t1,t2,displ1,displ2,gamma,height,half_duration_triangle

  logical :: triangle

  double precision, dimension(:), allocatable :: time,sem,sem_fil

! read file with number of lines in input
  open(unit=33,file='input_convolve_code.txt',status='old',action='read')
  read(33,*) nlines
  read(33,*) half_duration_triangle
  read(33,*) triangle
  close(33)

! allocate arrays
  allocate(time(nlines),sem(nlines),sem_fil(nlines))

! read the input seismogram
  do i = 1,nlines
    read(5,*) time(i),sem(i)
  enddo

! define a Gaussian with the right exponent to mimic a triangle of equivalent half duration
  alpha = SOURCE_DECAY_MIMIC_TRIANGLE/half_duration_triangle

! compute the time step
  dt = time(2) - time(1)

! number of integers for which the source wavelet is different from zero
  if(triangle) then
    N_j = ceiling(half_duration_triangle/dt)
  else
    N_j = ceiling(1.5d0*half_duration_triangle/dt)
  endif

  do i = 1,nlines

    sem_fil(i) = 0.d0

    do j = -N_j,N_j

      if(i > j .and. i-j <= nlines) then

      tau_j = dble(j)*dt

! convolve with a triangle
    if(triangle) then
       height = 1.d0 / half_duration_triangle
       if(abs(tau_j) > half_duration_triangle) then
         source = 0.d0
       else if (tau_j < 0.d0) then
         t1 = - N_j * dt
         displ1 = 0.d0
         t2 = 0.d0
         displ2 = height
         gamma = (tau_j - t1) / (t2 - t1)
         source= (1.d0 - gamma) * displ1 + gamma * displ2
       else
         t1 = 0.d0
         displ1 = height
         t2 = + N_j * dt
         displ2 = 0.d0
         gamma = (tau_j - t1) / (t2 - t1)
         source= (1.d0 - gamma) * displ1 + gamma * displ2
       endif

      else

! convolve with a Gaussian
        exponent = alpha**2 * tau_j**2
        if(exponent < 50.d0) then
          source = alpha*exp(-exponent)/sqrt(PI)
        else
          source = 0.d0
        endif

      endif

      sem_fil(i) = sem_fil(i) + sem(i-j)*source*dt

      endif

    enddo
  enddo

! compute number of samples to remove from end of seismograms
  number_remove = N_j + 1
  do i=1,nlines - number_remove
    write(*,*) sngl(time(i)),' ',sngl(sem_fil(i))
  enddo

  end program convolve_source_time_function

