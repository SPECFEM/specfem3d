
  program create_Hamming_window

  implicit none

  integer, parameter :: NSOURCES = 1000

  double precision, parameter :: xs_min = 0.020d0
  double precision, parameter :: xs_size = 0.010d0

  double precision, parameter :: factor_max = 1.d10

! pi
  double precision, parameter :: PI = 3.141592653589793d0

  integer :: isource

  double precision :: x,hamming

  do isource = 1,NSOURCES

! Hamming apodization window
! see e.g. http://docs.scipy.org/doc/numpy/reference/generated/numpy.hamming.html
! and http://www.mathworks.fr/fr/help/signal/ref/hamming.html
    x = dble(isource - 1) / dble(NSOURCES - 1)
    hamming = 0.54d0 - 0.46d0*cos(2*PI*x)

    write(*,*) '# source ',isource
    write(*,*) 'source_surf                     = .false.'
    write(*,*) 'xs                              = ',xs_min + x * xs_size
    write(*,*) 'zs                              = 0.075'
    write(*,*) 'source_type                     = 1'
    write(*,*) 'time_function_type              = 1'
    write(*,*) 'name_of_source_file             = YYYYYYYYYYYYYYYYYY'
    write(*,*) 'burst_band_width                = 0.'
    write(*,*) 'f0                              = 0.5d6'
    write(*,*) 'tshift                          = 0.d0'
    write(*,*) 'angleforce                      = 0.0'
    write(*,*) 'Mxx                             = 1.'
    write(*,*) 'Mzz                             = 1.'
    write(*,*) 'Mxz                             = 0.'
    write(*,*) 'factor                          = ',factor_max * hamming

  enddo

  end program create_Hamming_window

