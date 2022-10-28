!program test_deconv
!  implicit none
!  double complex , allocatable ::  RF(:)
!  double precision, allocatable :: UIN(:),WIN(:),gauss(:),time(:),RW(:)
!  double precision, allocatable :: U(:),W(:),P0(:),Radial(:),Vertical(:),P(:)
!  real(kind(0d0)), allocatable :: a(:)
!  complex(kind(0d0)), allocatable :: c(:)
!  double precision dt,f0,minderr,rms,t0,t1
!  double precision TSHIFT
!  integer nt,nt_rf,itnmax,i0,i1,i,nstep_convole_gauss
!  integer nfft
!  double precision exposant
!  integer exposant_int
!
!  f0=0.5d0
!  minderr=1.D-5
!  itnmax=200
!  TSHIFT=0.D0
!  nt=7700
!  dt=0.0125
!  t0=11.d0
!  t1=61.d0
!  allocate(Radial(nt),Vertical(nt),time(nt))


!  open(10,file='/home/monteill/Work/Specfem/Synth/DATA/data_sud_butt_0.001_0.95/S20.TS.HXX.semd')
!  read(10,*) time(1),Radial(1)
!  do i=2,nt
!     read(10,*)  time(i),Radial(i)
!     if (time(i-1) >= t0 .and. time(i) <= t0+2*dt ) i0=i
!     if (time(i-1) >= t1 .and. time(i) <= t1+2*dt ) i1=i
!  enddo
!  close(10)

!  open(10,file='/home/monteill/Work/Specfem/Synth/DATA/data_sud_butt_0.001_0.95/S20.TS.HXZ.semd')
!  do i=1,nt
!     read(10,*)  time(i),Vertical(i)
!  enddo
!  close(10)

!  nt_rf=i1-i0+1
  !allocate(RF(nt_rf),UIN(nt_rf),WIN(nt_rf),U(nt_rf),W(nt_rf),P0(nt_rf),RW(2*nt_rf),P(nt_rf))

!! test de la fft
!  exposant=(log(real(nt_rf,8))/log(2.d0))
!  exposant_int=ceiling(exposant)
!  nfft=2**exposant_int
!  write(*,*) 'nfft ', nfft
!  allocate(RF(2*nfft))
!  allocate(U(nfft))
!  U(:)=0.
!  write(*,*) i0,i1
!  U(1:nt_rf)=Radial(i0:i1)

!  RF(:)=dcmplx(0.d0)
!  RF(1:nt_rf)=dcmplx(U(:))
!  call fft(RF,nfft,1)
!  open(10,file='four1_test.txt')
!  do i=1,nfft
!     write(10,*) U(i),real(RF(i)),aimag(RF(i))
!  enddo



!!
!!$  nstep_convole_gauss = int(10.d0/(f0*dt))
!!$  write(*,*)  nstep_convole_gauss
!!$  allocate(gauss(nstep_convole_gauss))
!!$  call define_Gaussian(gauss,f0,dt,nstep_convole_gauss)
!!$
!!$  UIN(:)=Radial(i0:i1)
!!$  WIN(:)=vertical(i0:i1)
!!$
!!$  call makeRFdecon_la(RF,rms,UIN,WIN,dt,nt_rf,TSHIFT,f0,itnmax,minderr,gauss,U,W,P0,RW,P,nstep_convole_gauss)

!  open(10,file='RF.txt')
!  do i=1,nt_rf
!     write(10,*) RF(i)
!  enddo

!end program test_deconv

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
  ! je consideree un support compact (nstep_convole_function)
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

!-----

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


!############################################## DECONVOLUTION LIGORRIA & AMMON ###############
subroutine makeRFdecon_la(RF,rms,UIN,WIN,dt,nt,TSHIFT,f0,itnmax,minderr,gauss,U,W,P0,RW,P,nstep_convole_gauss)

  implicit none

  real(kind(0d0)) dt,TSHIFT,f0,rms,minderr
  integer nt,itnmax,nstep_convole_gauss
  real(kind(0.d0)) RF(nt),UIN(nt),WIN(nt),gauss(nstep_convole_gauss)
  real(kind(0.d0)) U(nt),W(nt),P0(nt),R(nt),RW(2*nt),P(nt)
  double precision powerU,tmp_var,sumsq_i,sumsq,derror,amp
  integer i,iter,ispike,find_max_abs


  write(*,*) f0/dt
  write(*,*) nstep_convole_gauss,f0,dt

  call convolve_src_function(dt,UIN,U,gauss,nt,nstep_convole_gauss)
  call convolve_src_function(dt,WIN,W,gauss,nt,nstep_convole_gauss)

  !open(10,file='gassian.txt')
  !do i=1,nt
  !   write(10,*) UIN(i),U(i)
  !enddo
  !close(10)


  P0(:)=0.d0
  R(:)=U(:)
  powerU=sum(U(:)*U(:))

  sumsq_i = 1.d0
  derror = 100.d0*powerU + minderr
  write(*,*) derror
  iter=0
  iter=itnmax-1
  do while (iter < itnmax .and. dabs(derror) > minderr)



    call convolve_function(dt,R,W,RW,nt)
    open(10,file='conv.txt')

    do i=1,2*nt
       write(10,*) RW(i)
    enddo
    close(10)

    open(10,file='vec_to_conv.txt')
    do i=1,nt
       write(10,*) R(i),W(i)
    enddo
    close(10)

    tmp_var=sum(W(:)*W(:))
    RW(:)=RW(:)/tmp_var

    ispike=find_max_abs(RW,2*nt)
    amp=RW(ispike)/dt
    P0(ispike)= P0(ispike) + amp
    write(*,*) 'iter :',iter, ' ispike :', ispike

    call convolve_src_function(dt,P0,P,gauss,nt,nstep_convole_gauss)
    call convolve_function(dt,P,W,RW,nt)

    R(:)=U(:) - RW(:)

    sumsq=sum(R(:)*R(:))/powerU
    derror=100.d0*(sumsq_i - sumsq)
    sumsq_i=sumsq

    iter=iter+1

  enddo

  call convolve_src_function(dt,P0,RF,gauss,nt,nstep_convole_gauss)



end subroutine makeRFdecon_la

subroutine convolve_function(dt,A,B,AB,n)
  implicit none
  integer n,i,j,j0,j1
  double precision dt,A(n),B(n),AB(2*n)

  AB(:)=0.d0

!!$  do i=1,n
!!$     do j=1,i
!!$        AB(i)=AB(i)+A(j)*B(n+1-j)*dt
!!$     enddo
!!$  enddo
!!$  do i=n+1,2*n
!!$     do j=i-n+1,n
!!$        AB(i)=AB(i)+A(j)*B(n+1-j)*dt
!!$        !if (i== 5000) write(*,*) j,n+1-j
!!$     enddo
!!$  enddo
  do i=1,2*n
     j0=max(1,i-n+1)
     j1=min(i,n)
     do j=j0,j1
        AB(i)=AB(i)+A(j)*B(i-j)
     enddo
  enddo

end subroutine convolve_function



function find_max_abs(X,n)
  implicit none
  integer find_max_abs,n,i
  double precision X(n), max_value

  max_value=0.d0

  do i=1,n
     if (max_value <= abs(X(i))) then
        find_max_abs=i
        max_value=X(i)
     endif
  enddo

end function

!-----------------------------------------

subroutine fft(dat, nn, isign)

  implicit none
  integer :: nn, isign
  real(kind(0d0)) :: dat(2*nn)
  real(kind(0d0)) :: wr, wi, wpr, wpi, wtemp, theta
  real(kind(0d0)) :: tempr, tempi
  integer :: n,i,j,m,mmax,istep

  n = 2*nn
  j = 1
  do i = 1,n,2
     if (j > i) then
        tempr = dat(j)
        tempi = dat(j+1)
        dat(j) = dat(i)
        dat(j+1) = dat(i+1)
        dat(i) = tempr
        dat(i+1) = tempi
     endif
     m = n/2

     do while ((m >= 2) .and. (j > m))
        j = j-m
        m = m/2
     enddo

     j = j+m

  enddo

  mmax = 2

  do while (n > mmax)
     istep = 2*mmax
     theta = 6.28318530717959d0/(isign*mmax)
     wpr = -2.d0*dsin(0.5d0*theta)**2
     wpi = dsin(theta)
     wr = 1.d0
     wi = 0.d0
     do m = 1, mmax, 2
        do i = m, n, istep
           j = i+mmax
           tempr = wr*dat(j) -wi*dat(j+1)
           tempi = wr*dat(j+1)+wi*dat(j)
           dat(j) = dat(i) - tempr
           dat(j+1) = dat(i+1) -tempi
           dat(i) = dat(i) + tempr
           dat(i+1) = dat(i+1) + tempi
        enddo
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
     enddo
     mmax = istep
  enddo
  return
end subroutine fft

!!$
!!$! subroutine FFT, Cooley-Tukey radix-2, Dif (decimation-in-frequency)
!!$! FFT program written by C. S. Burrus, Rice University, Sept. 1983.
!!$! Complex input data in arrays X (real part) and Y (imaginary part).
!!$      subroutine FFT(X,Y,N,M)
!!$      double precision X(1),Y(1)
!!$! Main FFT loops
!!$      N2=N
!!$      DO 10 K=1,M
!!$         N1=N2
!!$         N2=N2/2
!!$         E=6.28318531/N1
!!$         A=0.
!!$         DO 20 J=1,N2
!!$            C=COS(A)
!!$            S=SIN(A)
!!$            A=E*J
!!$            DO 30 I=J,N,N1
!!$               L=I+N2
!!$               XT=X(I)-X(L)
!!$               X(I)=X(I)+X(L)
!!$               YT=Y(I)-Y(L)
!!$               Y(I)=Y(I)+Y(L)
!!$               X(L)=C*XT+S*YT
!!$               Y(L)=C*YT-S*XT
!!$            enddo
!!$         enddo
!!$      enddo
!!$      ! Digit reverse counter
!!$100   J=1
!!$      N1=N-1
!!$      DO 104 I=1,N1
!!$         if (I >= J) goto 101
!!$         XT=X(J)
!!$         X(J)=X(I)
!!$         X(I)=XT
!!$         XT=Y(J)
!!$          Y(J)=Y(I)
!!$         Y(I)=XT
!!$  101    K=N/2
!!$102      if (K >= J) goto 103
!!$         J=J-K
!!$         K=K/2
!!$         goto 102
!!$103      J=J+K
!!$      enddo
!!$      return
!!$    end subroutine FFT
!  four1 -- this is the four1 routine from numerical recipes
      subroutine FOUR1(DATA,NN,ISIGN)
      double precision  WR,WI,WPR,WPI,WTEMP,THETA
      double complex DATA(2*NN)
      N=2*NN
      J=1
      DO  I=1,N,2
        if (J > I) then
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        endif
        M=N/2
1       if ((M >= 2) .and. (J > M)) then
          J=J-M
          M=M/2
        goto 1
        endif
        J=J+M
     enddo
      MMAX=2
2     if (N > MMAX) then
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO  M=1,MMAX,2
          DO  I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
         enddo
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
       enddo
        MMAX=ISTEP
      goto 2
      endif
      return
    end subroutine FOUR1

!
! -------------------------------------------------
!
