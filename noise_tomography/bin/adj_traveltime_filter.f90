program adj_traveltime

implicit none

integer, parameter :: nstep = 2999
double precision, parameter :: dt = 0.05d0
double precision, parameter :: C_crit = 0.5d0

double precision, parameter :: length_time_window = 70d0
integer, parameter :: taper_type = 1
integer, parameter :: i_skip_window =1
double precision, parameter :: factor = 5d-2


integer, parameter :: nrec = 1
double precision, parameter :: pi = 3.141592653589793d0

double precision :: data_origin(nstep,nrec),syn_origin(nstep,nrec),adj(nstep,nrec),adj_density(nstep,nrec)
double precision :: data_picked(nstep,nrec),syn_picked(nstep,nrec),data_filtered(nstep,nrec),syn_filtered(nstep,nrec)
double precision :: data_temp(nstep),syn_temp(nstep),data_direct(nstep),syn_direct(nstep)
double precision :: trace_data_max(nrec), trace_syn_max(nrec), t(nstep),data_max,syn_max
integer :: flag(nrec), irec, itime, length_window, i_start_window, i_end_window
character(len=256) :: station_name,file_data,file_syn,file_data_direct,file_syn_direct
character(len=256) :: file_adj,file_adj_density, file_misfit, file_adj_BXX,file_adj_BXZ, file_adj_zeros
double precision :: c(nstep), data_trace(nstep), syn_trace(nstep)
integer :: I(nstep)
double precision :: Norm_data_crit,Norm_syn_crit,AMP,Norm_data,Norm_syn,Norm,Norm_adj_temp
double precision :: c_current,c_final(100), w(nstep), misfit_traveltime, traveltime_delay
integer :: I_current, n_current, index_current, l_current, l, lag
integer :: I_final(100), index_final(100)

double precision, allocatable :: taper(:),corre(:), adj_temp(:),data_trace_temp(:),syn_trace_temp(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer,parameter :: filter_flag = 1
integer :: ifreq, nfreq
real :: F1,F2,D(8),G,DELT
!real  freq_low(2),freq_high(2)
!data  freq_low  / 1.0d-4 , 1.0d-4/
!data  freq_high / 5.0d-1 , 5.0d-2/
real  freq_low(1),freq_high(1)
data  freq_low  / 0.01d0 /
data  freq_high / 0.2d0 /
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nfreq=size(freq_low)

file_misfit = '../OUTPUT_FILES/SEM/misfit_traveltime_delay'
open(unit=1111,file=trim(file_misfit),status='unknown')

DELT=dt*1.0d3
adj=0.0d0
misfit_traveltime = 0.0d0
!!!! loading data and synthetics !!!!
do irec = 1,nrec
   file_data         = '../OUTPUT_FILES/SEM/X2.DB.BXZ.semd'
   open(unit=1001,file=trim(file_data),status='old',action='read')
   do itime = 1,nstep
           read(1001,*) t(itime),data_origin(nstep-itime+1,irec)  ! the reversed seismogram involves $C^\alpha\beta(t)=C^\beta\alpha(-t)$
   end do
   close(1001)

   data_temp(1)=0.0
   data_temp(nstep)=0.0
   do itime = 2,nstep-1
      data_temp(itime)=( data_origin(itime+1,irec) - data_origin(itime-1,irec) )/(2*dt)                  ! displacement --> particle  velocity (elastic  ".semd")
   end do
   data_origin(:,irec)=data_temp
   trace_data_max(irec)=maxval(abs(data_temp))
end do

data_max=maxval(trace_data_max)

!!!! taper !!!!
length_window = floor(length_time_window/dt)+1
allocate(taper(length_window)); taper(:)=1.0d0
allocate(adj_temp(length_window)); adj_temp(:)=0.0d0
allocate(corre(2*length_window-1)); corre(:)=0.0d0
allocate(data_trace_temp(length_window)); data_trace_temp(:)=0.0d0
allocate(syn_trace_temp(length_window)); syn_trace_temp(:)=0.0d0
if (taper_type .eq. 1) then  ! cosine taper, otherwise using a constant (1.0) instead
   do l=1,length_window
      taper(l) = (1.0-cos(pi*2.0*(l-1)/(length_window-1)))/2.0
   end do
end if

!!!! computing adj sources !!!!
                  if (filter_flag .eq. 0)   then
                      nfreq=1
                  end if
do irec = 1,nrec
   do ifreq = 1,nfreq
     ! bandpass filter
     data_filtered(:,irec)=data_origin(:,irec)
      syn_filtered(:,irec)= syn_origin(:,irec)
     if (filter_flag .eq. 1) then
         print *, 'filtering'
         write(*,"('frequency band ',i2,' : f_min = ',f10.5,' Hz , f_max = ',f10.5, ' Hz')") ifreq, freq_low(ifreq), freq_high(ifreq)
     ! THIS SECTION CALCULATES THE FILTER AND MUST BE CALLED BEFORE
     ! FILTER IS CALLED
         F1=freq_low(ifreq)
         F2=freq_high(ifreq)
         call BNDPAS(F1,F2,DELT,D,G,nstep)
     !    F1 = LOW FREQUENCY CUTOFF (6 DB DOWN)
     !    F2 = HIGH FREQUENCY CUTOFF (6 DB DOWN)
     !    DELT = SAMPLE INTERVAL IN MILLISECONDS
     !    D = WILL CONTAIN 8 Z DOMAIN COEFICIENTS OF RECURSIVE FILTER
     !    G = WILL CONTAIN THE GAIN OF THE FILTER,
         call FILTER(data_filtered(:,irec),nstep,D,G,2)
         call FILTER(syn_filtered(:,irec),nstep,D,G,2)
     !     X = DATA VECTOR OF LENGTH N CONTAINING DATA TO BE FILTERED
     !     D = FILTER COEFFICIENTS CALCULATED BY BNDPAS
     !     G = FILTER GAIN
     !     IG = 1  one pass
     !     ig = 2  two passes
     end if

     c(:)=0.0
     I(:)=(length_window-1)
     data_trace=data_filtered(:,irec)
      syn_trace= syn_filtered(:,irec)
     !!!! cross-correlation !!!!
     do i_start_window = (nstep+1)/2, (nstep+1)/2                                  ! the positive branch
     !do i_start_window = (nstep+1)/2 - length_window, (nstep+1)/2 - length_window  ! the negative branch
        i_end_window = i_start_window + length_window - 1

        data_trace_temp=data_trace(i_start_window:i_end_window)

        Norm_data_crit=sqrt(DOT_PRODUCT(data_trace,data_trace))
        AMP = Norm_data_crit

        Norm_data=sqrt(DOT_PRODUCT( data_trace_temp,data_trace_temp ))
        Norm = Norm_data*Norm_data

        w(:)=0.0
        w(i_start_window : i_end_window)=taper
     end do

     data_trace=data_trace*w
     data_filtered(:,irec)=data_trace
     !!!! normal adjoint sources !!!!
           adj_temp=data_filtered(i_start_window : i_end_window,irec)
           Norm_adj_temp = - DOT_PRODUCT(adj_temp,adj_temp)*dt  ! minus sign comes from integration by part !
           adj_temp=adj_temp/Norm_adj_temp                                            ! ray density map/banana-doughnut kernel (DeltaT=+1)
           adj(i_start_window : i_end_window ,irec )=adj_temp
           traveltime_delay=dt
           misfit_traveltime = misfit_traveltime + traveltime_delay * traveltime_delay / 2.0
           write(1111,*) traveltime_delay

   end do  !do ifreq=1,nfreq

   !!!! output !!!!
   file_adj_BXZ      = '../OUTPUT_FILES/SEM/adj_sources_contribution1'
   open(unit=1002,file=trim(file_adj_BXZ),status='unknown')
   do itime = 1,nstep
      write(1002,*) t(itime), adj(nstep-itime+1,irec)
   end do
   close(1002)
   file_adj_zeros      = '../OUTPUT_FILES/SEM/adj_sources_contribution2'
   open(unit=1002,file=trim(file_adj_zeros),status='unknown')
   do itime = 1,nstep
      write(1002,*) t(itime), adj(itime,irec)
   end do
   close(1002)

enddo  !do irec=1,nrec


close(1111)

file_misfit = '../OUTPUT_FILES/SEM/misfit_traveltime'
open(unit=1001,file=trim(file_misfit),status='unknown')
write(1001,*) misfit_traveltime, traveltime_delay
close(1001)

end program adj_traveltime










!!!!!!! subroutine for bandpass filter, from some website !!!!!!!!!
! http://www-lgit.obs.ujf-grenoble.fr/users/jrevilla/seiscomp/patch/pack/plugins/seisan/LIB/bndpas.for

SUBROUTINE BNDPAS(F1,F2,DELT,D,G,N)
! RECURSIVE BUTTERWORTH BAND PASS FILTER (KANASEWICH, TIME SERIES
! ANALYSIS IN GEOPHYSICS, UNIVERSITY OF ALBERTA PRESS, 1975; SHANKS,
! JOHN L, RECURSION FILTERS FOR DIGITAL PROCESSING, GEOPHYSICS, V32,
! FILTER.  THE FILTER WILL HAVE 8 POLES IN THE S PLANE AND IS
! APPLIED IN FORWARD AND REVERSE DIRECTIONS SO AS TO HAVE ZERO
! PHASE SHIFT.  THE GAIN AT THE TWO FREQUENCIES SPECIFIED AS
! CUTOFF FREQUENCIES WILL BE -6DB AND THE ROLLOFF WILL BE ABOUT
! THE FILTER TO PREVENT ALIASING PROBLEMS.
    COMPLEX P(4),S(8),Z1,Z2
    real D(8),XC(3),XD(3),XE(3)
    double precision :: X(N)
    DATA ISW/0/,TWOPI/6.2831853/
! THIS SECTION CALCULATES THE FILTER AND MUST BE CALLED BEFORE
! FILTER IS CALLED

!    F1 = LOW FREQUENCY CUTOFF (6 DB DOWN)
!    F2 = HIGH FREQUENCY CUTOFF (6 DB DOWN)
!    DELT = SAMPLE INTERVAL IN MILLISECONDS
!    D = WILL CONTAIN 8 Z DOMAIN COEFICIENTS OF RECURSIVE FILTER
!    G = WILL CONTAIN THE GAIN OF THE FILTER,

      DT=DELT/1000.0
      TDT=2.0/DT
      FDT=4.0/DT
      ISW=1
      P(1)=CMPLX(-.3826834,.9238795)
      P(2)=CMPLX(-.3826834,-.9238795)
      P(3)=CMPLX(-.9238795,.3826834)
      P(4)=CMPLX(-.9238795,-.3826834)
      W1=TWOPI*F1
      W2=TWOPI*F2
      W1=TDT*TAN(W1/TDT)
      W2=TDT*TAN(W2/TDT)
      HWID=(W2-W1)/2.0
      WW=W1*W2
      DO 19 I=1,4
      Z1=P(I)*HWID
      Z2=Z1*Z1-WW
      Z2=CSQRT(Z2)
      S(I)=Z1+Z2
   19 S(I+4)=Z1-Z2
      G=.5/HWID
      G=G*G
      G=G*G
      DO 29 I=1,7,2
      B=-2.0*REAL(S(I))
      Z1=S(I)*S(I+1)
      C=REAL(Z1)
      A=TDT+B+C/TDT
      G=G*A
      D(I)=(C*DT-FDT)/A
   29 D(I+1)=(A-2.0*B)/A
      G=G*G
    5 FORMAT ('-FILTER GAIN IS ', 9E12.6)
      RETURN

      ENTRY FILTER(X,N,D,G,IG)

!print*, maxval(X), ISW
!     X = DATA VECTOR OF LENGTH N CONTAINING DATA TO BE FILTERED
!     D = FILTER COEFFICIENTS CALCULATED BY BNDPAS
!     G = FILTER GAIN
!     IG = 1  one pass
!     ig = 2  two passes

      IF (ISW.EQ.1) GO TO 31
      WRITE (6,6)
    6 FORMAT ('1BNDPAS MUST BE CALLED BEFORE FILTER')
      return

!     APPLY FILTER IN FORWARD DIRECTION

   31 XM2=X(1)
      XM1=X(2)
      XM=X(3)
      XC(1)=XM2
      XC(2)=XM1-D(1)*XC(1)
      XC(3)=XM-XM2-D(1)*XC(2)-D(2)*XC(1)
      XD(1)=XC(1)
      XD(2)=XC(2)-D(3)*XD(1)
      XD(3)=XC(3)-XC(1)-D(3)*XD(2)-D(4)*XD(1)
      XE(1)=XD(1)
      XE(2)=XD(2)-D(5)*XE(1)
      XE(3)=XD(3)-XD(1)-D(5)*XE(2)-D(6)*XE(1)
      X(1)=XE(1)
      X(2)=XE(2)-D(7)*X(1)
      X(3)=XE(3)-XE(1)-D(7)*X(2)-D(8)*X(1)
      DO 39 I=4,N
      XM2=XM1
      XM1=XM
      XM=X(I)
      K=I-((I-1)/3)*3
      GO TO (34,35,36),K
   34 M=1
      M1=3
      M2=2
      GO TO 37
   35 M=2
      M1=1
      M2=3
      GO TO 37
   36 M=3
      M1=2
      M2=1
   37 XC(M)=XM-XM2-D(1)*XC(M1)-D(2)*XC(M2)
      XD(M)=XC(M)-XC(M2)-D(3)*XD(M1)-D(4)*XD(M2)
      XE(M)=XD(M)-XD(M2)-D(5)*XE(M1)-D(6)*XE(M2)
   39 X(I)=XE(M)-XE(M2)-D(7)*X(I-1)-D(8)*X(I-2)
!
!
      if(ig.eq.1) goto 3333
      XM2=X(N)
      XM1=X(N-1)
      XM=X(N-2)
      XC(1)=XM2
      XC(2)=XM1-D(1)*XC(1)
      XC(3)=XM-XM2-D(1)*XC(2)-D(2)*XC(1)
      XD(1)=XC(1)
      XD(2)=XC(2)-D(3)*XD(1)
      XD(3)=XC(3)-XC(1)-D(3)*XD(2)-D(4)*XD(1)
      XE(1)=XD(1)
      XE(2)=XD(2)-D(5)*XE(1)
      XE(3)=XD(3)-XD(1)-D(5)*XE(2)-D(6)*XE(1)
      X(N)=XE(1)
      X(N-1)=XE(2)-D(7)*X(1)
      X(N-2)=XE(3)-XE(1)-D(7)*X(2)-D(8)*X(1)
      DO 49 I=4,N
      XM2=XM1
      XM1=XM
      J=N-I+1
      XM=X(J)
      K=I-((I-1)/3)*3
      GO TO (44,45,46),K
   44 M=1
      M1=3
      M2=2
      GO TO 47
   45 M=2
      M1=1
      M2=3
      GO TO 47
   46 M=3
      M1=2
      M2=1
   47 XC(M)=XM-XM2-D(1)*XC(M1)-D(2)*XC(M2)
      XD(M)=XC(M)-XC(M2)-D(3)*XD(M1)-D(4)*XD(M2)
      XE(M)=XD(M)-XD(M2)-D(5)*XE(M1)-D(6)*XE(M2)
   49 X(J)=XE(M)-XE(M2)-D(7)*X(J+1)-D(8)*X(J+2)
 3333 continue
      if (ig.eq.1) then
        gg=sqrt(g)   ! if only pass once, modify gain
      else
        gg=g
      endif
      DO 59 I=1,N
   59 X(I)=X(I)/gg
      RETURN
END

