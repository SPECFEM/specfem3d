!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, April 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
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

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0

  integer :: i,j,N_j,nlines !!!!!!!!! ,number_remove

  double precision :: dt,tau_j,source,exponent

  double precision, dimension(:), allocatable :: time,sem,sem_fil

! other parameters for the source
  double precision :: a,f0,sigma_source

! read file with number of lines in input
  open(unit=33,file='input_convolve_code.txt',status='old',action='read')
  read(33,*) nlines
  close(33)

! allocate arrays
  allocate(time(nlines),sem(nlines),sem_fil(nlines))

! read the input seismogram
  do i = 1,nlines
! ignore 1st and 3rd components in Li Zhao's file format
    read(5,*) time(i),sem(i)
  enddo

! compute the time step
  dt = time(2) - time(1)

!! DK DK define the Gaussian source
  sigma_source = 1.8d0 !! 2.d0 !!! 1.5d0 !!!!!  2.5d0 !!!! 5.d0 !!!!!!!!! 30.d0 * 0.8d0 * DT
  f0 = 1.d0 / (PI * sigma_source * sqrt(2.d0))

! number of integers for which the source wavelet is different from zero
!!!!  N_j = 20.d0 / dt
  N_j = 100.d0 / dt
!!!  N_j = nlines + 10
!!!!!!!!   N_j = ceiling(1.5d0*half_duration_triangle/dt)

  do i = 1,nlines

    sem_fil(i) = 0.d0

    do j = -N_j,N_j

      if(i > j .and. i-j <= nlines) then

      tau_j = dble(j)*dt

!! DK DK put a small Gaussian for the source to mimic a Dirac
!!!!!!!!!!!!!!!!!      time_source = (it-1)*DT
      a = pi*pi*f0*f0
!!!!!!!!!!!!!      stf = exp(-a*(time_source-t0)**2)

! convolve with a Gaussian
!!!!!!!!!!!!!!!        exponent = alpha**2 * tau_j**2
        exponent = a * tau_j**2
        if(exponent < 50.d0) then
          source = exp(-exponent)
        else
          source = 0.d0
        endif

      sem_fil(i) = sem_fil(i) + sem(i-j)*source*dt

      endif

    enddo
  enddo

! compute number of samples to remove from end of seismograms
!!!!!!!!!!!!!!!!  number_remove = N_j + 1
  do i=1,nlines !!!!!!!!!!!!!!!! - number_remove
    write(*,*) sngl(time(i)),sngl(sem_fil(i)) !!!!!!!!!!! * 4000000000.d0)
  enddo

  end program convolve_source_time_function

