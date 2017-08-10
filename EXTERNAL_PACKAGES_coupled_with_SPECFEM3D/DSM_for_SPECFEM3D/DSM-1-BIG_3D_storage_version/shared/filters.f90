!###################################### PREMIER FILTRE #################################################
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

  real(kind(0d0)), dimension(1)::x,y
  real(kind(0d0)), dimension (10) ::  a, b1, b2
  real(kind(0d0)) :: dt,f1,f2
  integer :: iunit, npoles,norder,irek,n,lx
  !real(kind(0d0)) :: x(n),y(n)

   iunit = 3

   if (norder /= 0) then
      npoles=iabs(norder)
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

!---------------------------------------------------------------


subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
  ! performs recursive filtering of data in array x of length ndat
  ! filtered output in y
  ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
  ! npoles is the number of poles
  ! iflag=0: forward filtering only
  ! iflag /= 0: forward and backward filtering

  implicit none

  real(kind(0d0)), dimension(10) :: z,z1,z2 ,a,b1,b2
  real(kind(0d0)) ::  x1,x2
  integer :: ndat, npoles, iflag, n,i
  real(kind(0d0)) :: x(ndat), y(ndat)

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



!---------------------------------------------------------------


subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
  !determines filtercoefficients for recursive bandpassfilter

  real(kind(0d0)),dimension(10) :: a,b1,b2
  complex(kind(0d0)) :: s(20), t1,t2,p
  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)) :: f1,f2,dt,d2,w0,w1,w2,ssum, sprod,fact1,fact2,fact3
  integer :: i,npol2,n,npoles


  if (npoles > 10) then
     stop ' npoles greater than 10: STOP '
  endif

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
  return
end subroutine bpcoeff




!####################### FIN PREMEIER ######################################################################################

!############################## DEUXIEME FILTRE ##############################
subroutine convolve_src_function(dt,y_input,y_output,convole_function,nsamples,nstep_convole_function)
  implicit none
  integer nsamples,nstep_convole_function
  real(kind(0d0)) dt,y_input(nsamples),y_output(nsamples),time(nsamples),convole_function(nstep_convole_function)
  integer it,j,jmin,jmax

  y_output(:)=0.d0
  do it=1,nsamples
     jmin=max(1,ceiling(it-0.5*(nstep_convole_function-1)))
     jmax=min(nsamples,floor(it+0.5*(nstep_convole_function-1)))
     do j=jmin,jmax
        y_output(it) = y_output(it) + y_input(j)*convole_function(1+j-jmin)*dt
     enddo
  enddo

end subroutine convolve_src_function

subroutine define_Gaussian(gauss,f0,dt,n)
 implicit none
 integer i,n
 real(kind(0d0)) dt,f0,gauss(n),i0,const
 const=f0/sqrt(3.141592653589793)
 i0=n/2
 do i=1,n
   gauss(i)=const*exp(- ((i - i0)*dt*f0)**2)
 enddo
end subroutine define_Gaussian
!####################################### FIN DEUXIEME ####################
