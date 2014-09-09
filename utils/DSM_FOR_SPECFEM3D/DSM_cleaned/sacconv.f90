
  subroutine convspc ( max_nstation,maxnfreq,n_station, &
           time_series_length,n_frequency,omega_imag, &
           station_displacement, &
           source_depth,source_lat,source_lon,station_lat,station_lon, &
           work_spc,work_spc_real,work_time,sac_file,ITYPE_SEISMOGRAMS )

  implicit none

! variables for input/output
  integer max_nstation,maxnfreq,n_frequency,n_station,ITYPE_SEISMOGRAMS
  real(kind=8) time_series_length,omega_imag
  real(kind=8) source_depth,source_lat,source_lon
  real(kind=8) station_lat(max_nstation),station_lon(max_nstation)
  complex(kind=8) station_displacement(3,max_nstation,maxnfreq)
  real work_time(32*maxnfreq)
  character(len=80) sac_file(max_nstation)

  complex(kind=8) work_spc(16*maxnfreq)
!! DK DK added this to avoid a incompatible type warning when calling routine the Fourier transform routine later
!! DK DK the two arrays work_spc and work_spc_real are equivalenced in the calling program,
!! DK DK and thus they correspond to the same memory block
  real(kind=8) work_spc_real(2*16*maxnfreq)

! other variables
  integer i_station,icomp,i_frequency,ismooth,n1,m1,nn,isignval,i
  real(kind=8) t
  character(len=4) cext(3)
  character(len=80) outfile

  real(kind=8), parameter :: pi=3.1415926535897932d0

  cext(1) = '.bhz'
  cext(2) = '.bhr'
  cext(3) = '.bht'

  ismooth = 8

  do i_station=1,n_station

  do icomp=1,3

    work_spc(1) = cmplx(0.d0)

! the code computes velocity Green functions for Heaviside function source time history
! (or displacement Green functions for delta function source time history).
! To obtain displacement seismograms, we need to integrate the obtained seismograms.
    if(ITYPE_SEISMOGRAMS == 2) then   !   compute velocity
      do i_frequency=1,n_frequency
        work_spc(i_frequency+1) = station_displacement(icomp,i_station,i_frequency)
      enddo

    else if(ITYPE_SEISMOGRAMS == 1) then   !   compute displacement
      do i_frequency=1,n_frequency
        ! divide the spectrum by i\omega
        work_spc(i_frequency+1) = station_displacement(icomp,i_station,i_frequency) &
              / cmplx( 0.d0, 2.d0 * pi * dble(i_frequency) / time_series_length)
      enddo

    else
      stop 'error: incorrect value of ITYPE_SEISMOGRAMS, should be 1 or 2'
    endif

    do i_frequency=n_frequency+2,ismooth*n_frequency+1
      work_spc(i_frequency) = cmplx(0.d0)
    enddo

    do i_frequency=1,ismooth*n_frequency-1
      n1 = ismooth*n_frequency + i_frequency + 1
      m1 = ismooth*n_frequency - i_frequency + 1
      work_spc(n1) = conjg( work_spc(m1) )
    enddo

    nn = 2 * ismooth*n_frequency

    isignval = 1
!! DK DK the two arrays work_spc and work_spc_real are equivalenced in the calling program,
!! DK DK and thus they correspond to the same memory block
!! DK DK I thus replace one with the other here to avoid a warning by the compiler about incompatible types in the subroutine
!! DK DK    call fast_fourier_transform_real( work_spc,nn,isignval )
    call fast_fourier_transform_real( work_spc_real,nn,isignval )

    do i=-nn+1,nn
      t = time_series_length * dble(i-1) / dble(nn)
      work_time(i+nn) = real( dble( work_spc( mod(i-1+nn,nn)+1 ) ) * dexp( omega_imag * t ) )

! amplitude correction
! --- for FFT
      work_time(i+nn) = work_time(i+nn) / real( time_series_length )

! --- convert kilometers to meters
      work_time(i+nn) = work_time(i+nn) * 1.e3

    enddo
    outfile = sac_file(i_station)(1:len_trim(sac_file(i_station))) //cext(icomp)
    call sacconv( 2*time_series_length,2*nn, -time_series_length,work_time, &
                          source_depth,source_lat,source_lon, &
                          station_lat(i_station),station_lon(i_station),icomp,outfile )
    enddo
  enddo

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sacconv( tlen,np,b0,y,d0,theta0,phi0,theta,phi,icomp,output )

! conversion of the seismograms to the Seismic Analysis Code (SAC) binary format

  implicit none

  real(kind=8), parameter :: pi=3.1415926535897932d0

! header variables (I)
  integer, parameter :: mfhdr=70, mnhdr=15, mihdr=20, mlhdr=5, mkhdr=24
  real fhdr(mfhdr)
  integer nhdr(mnhdr),ihdr(mihdr),lhdr(mlhdr)
  character(len=8) khdr(mkhdr)

! header variables (II)
  real delta,depmin,depmax,depmn,depmx,scaleval,odelta
  real b,begin,e,ennd,o,origin,a,arrivl,fmtval
  real t0,t1,t2,t3,t4,t5,t6,t7,t8,t9
  real time0,time1,time2,time3,time4
  real time5,time6,time7,time8,time9
  real f,fini
  real resp0,resp1,resp2,resp3,resp4
  real resp5,resp6,resp7,resp8,resp9
  real stla,stlo,stel,stdp,evla,evlo,evel,evdp,fhdr40
  real user0,user1,user2,user3,user4
  real user5,user6,user7,user8,user9
  real dist,az,baz,gcarc,sb,sdelta,depmen,fmean,cmpaz,cmpinc
  real xminimum,xmaximum,yminimum,ymaximum
  real fhdr64,fhdr65,fhdr66,fhdr67,fhdr68,fhdr69,fhdr70
  integer nzdttm(6),nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec
  integer nvhdr,ninf,nhst,npts,nsnpts,nsn,nxsize,nysize,nhdr15
  integer iftype,idep,iztype,ihdr4,iinst
  integer istreg,ievreg,ievtyp,iqual,isynth
  integer ihdr11,ihdr12,ihdr13,ihdr14,ihdr15
  integer ihdr16,ihdr17,ihdr18,ihdr19,ihdr20
  integer leven,lpspol,lovrok,lcalda,lhdr5

  character(len=8) kstnm
  character(len=16) kevnm
  character(len=8) khole,ko,ka,kt0,kt1,kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9
  character(len=8) kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,kinst

! equivalence declaration
  equivalence (delta,fhdr(01))
  equivalence (depmin,depmn,fhdr(02))
  equivalence (depmax,depmx,fhdr(03))
  equivalence (scaleval,fhdr(04))
  equivalence (odelta,fhdr(05))
  equivalence (b,begin,fhdr(06))
  equivalence (e,ennd,fhdr(07))
  equivalence (o,origin,fhdr(08))
  equivalence (a,arrivl,fhdr(09))
  equivalence (fmtval,fhdr(10))
  equivalence (t0,time0,fhdr(11))
  equivalence (t1,time1,fhdr(12))
  equivalence (t2,time2,fhdr(13))
  equivalence (t3,time3,fhdr(14))
  equivalence (t4,time4,fhdr(15))
  equivalence (t5,time5,fhdr(16))
  equivalence (t6,time6,fhdr(17))
  equivalence (t7,time7,fhdr(18))
  equivalence (t8,time8,fhdr(19))
  equivalence (t9,time9,fhdr(20))
  equivalence (f,fini,fhdr(21))
  equivalence (resp0,fhdr(22))
  equivalence (resp1,fhdr(23))
  equivalence (resp2,fhdr(24))
  equivalence (resp3,fhdr(25))
  equivalence (resp4,fhdr(26))
  equivalence (resp5,fhdr(27))
  equivalence (resp6,fhdr(28))
  equivalence (resp7,fhdr(29))
  equivalence (resp8,fhdr(30))
  equivalence (resp9,fhdr(31))
  equivalence (stla,fhdr(32))
  equivalence (stlo,fhdr(33))
  equivalence (stel,fhdr(34))
  equivalence (stdp,fhdr(35))
  equivalence (evla,fhdr(36))
  equivalence (evlo,fhdr(37))
  equivalence (evel,fhdr(38))
  equivalence (evdp,fhdr(39))
  equivalence (fhdr40,fhdr(40))
  equivalence (user0,fhdr(41))
  equivalence (user1,fhdr(42))
  equivalence (user2,fhdr(43))
  equivalence (user3,fhdr(44))
  equivalence (user4,fhdr(45))
  equivalence (user5,fhdr(46))
  equivalence (user6,fhdr(47))
  equivalence (user7,fhdr(48))
  equivalence (user8,fhdr(49))
  equivalence (user9,fhdr(50))
  equivalence (dist,fhdr(51))
  equivalence (az,fhdr(52))
  equivalence (baz,fhdr(53))
  equivalence (gcarc,fhdr(54))
  equivalence (sb,fhdr(55))
  equivalence (sdelta,fhdr(56))
  equivalence (depmen,fmean,fhdr(57))
  equivalence (cmpaz,fhdr(58))
  equivalence (cmpinc,fhdr(59))
  equivalence (xminimum,fhdr(60))
  equivalence (xmaximum,fhdr(61))
  equivalence (yminimum,fhdr(62))
  equivalence (ymaximum,fhdr(63))
  equivalence (fhdr64,fhdr(64))
  equivalence (fhdr65,fhdr(65))
  equivalence (fhdr66,fhdr(66))
  equivalence (fhdr67,fhdr(67))
  equivalence (fhdr68,fhdr(68))
  equivalence (fhdr69,fhdr(69))
  equivalence (fhdr70,fhdr(70))
  equivalence (nzyear,nzdttm,nhdr(1))
  equivalence (nzjday,nhdr(2))
  equivalence (nzhour,nhdr(3))
  equivalence (nzmin,nhdr(4))
  equivalence (nzsec,nhdr(5))
  equivalence (nzmsec,nhdr(6))
  equivalence (nvhdr,nhdr(7))
  equivalence (ninf,nhdr(8))
  equivalence (nhst,nhdr(9))
  equivalence (npts,nhdr(10))
  equivalence (nsnpts,nhdr(11))
  equivalence (nsn,nhdr(12))
  equivalence (nxsize,nhdr(13))
  equivalence (nysize,nhdr(14))
  equivalence (nhdr15,nhdr(15))
  equivalence (iftype,ihdr(1))
  equivalence (idep,ihdr(2))
  equivalence (iztype,ihdr(3))
  equivalence (ihdr4,ihdr(4))
  equivalence (iinst,ihdr(5))
  equivalence (istreg,ihdr(6))
  equivalence (ievreg,ihdr(7))
  equivalence (ievtyp,ihdr(8))
  equivalence (iqual,ihdr(9))
  equivalence (isynth,ihdr(10))
  equivalence (ihdr11,ihdr(11))
  equivalence (ihdr12,ihdr(12))
  equivalence (ihdr13,ihdr(13))
  equivalence (ihdr14,ihdr(14))
  equivalence (ihdr15,ihdr(15))
  equivalence (ihdr16,ihdr(16))
  equivalence (ihdr17,ihdr(17))
  equivalence (ihdr18,ihdr(18))
  equivalence (ihdr19,ihdr(19))
  equivalence (ihdr20,ihdr(20))
  equivalence (leven,lhdr(1))
  equivalence (lpspol,lhdr(2))
  equivalence (lovrok,lhdr(3))
  equivalence (lcalda,lhdr(4))
  equivalence (lhdr5,lhdr(5))
  equivalence (kstnm,khdr(1))
  equivalence (kevnm,khdr(2))
  equivalence (khole,khdr(4))
  equivalence (ko,khdr(5))
  equivalence (ka,khdr(6))
  equivalence (kt0,khdr(7))
  equivalence (kt1,khdr(8))
  equivalence (kt2,khdr(9))
  equivalence (kt3,khdr(10))
  equivalence (kt4,khdr(11))
  equivalence (kt5,khdr(12))
  equivalence (kt6,khdr(13))
  equivalence (kt7,khdr(14))
  equivalence (kt8,khdr(15))
  equivalence (kt9,khdr(16))
  equivalence (kf,khdr(17))
  equivalence (kuser0,khdr(18))
  equivalence (kuser1,khdr(19))
  equivalence (kuser2,khdr(20))
  equivalence (kcmpnm,khdr(21))
  equivalence (knetwk,khdr(22))
  equivalence (kdatrd,khdr(23))
  equivalence (kinst,khdr(24))

  integer np,icomp,i,ier
  real(kind=8) tlen,b0,d0
  real y(*)
  real(kind=8) theta0,phi0,theta,phi

! variables for files
  character(len=80) output

  integer, external :: filenchk,sacbw

! initialize the SAC header
  do i=1,mfhdr
    fhdr(i) = -12345.0
  enddo

  do i=1,mnhdr
    nhdr(i) = -12345
  enddo

  do i=1,mihdr
    ihdr(i) = -12345
  enddo

  do i=1,mlhdr
    lhdr(i) = 0
  enddo

  do i=1,mkhdr
    khdr(i) = '-12345  '
  enddo

  begin = real(b0)
  npts = np
  delta = tlen / real(npts)
  stla = real( theta )
  stlo = real( phi )
  evdp = real( d0 )
  evla = real( theta0 )
  evlo = real( phi0 )

! generate SAC-binary file for evenly sampled data

!     icomp
!       3: Transverse Component
!       otherwise: not supported by this version

  depmax = y(1)
  depmin = y(2)
  depmen = 0.d0

  do i=1,npts
    if ( y(i)>depmax ) depmax = y(i)
    if ( y(i)<depmin ) depmin = y(i)
    depmen = depmen + y(i)
  enddo

  depmen = depmen / real(npts)

  ennd = begin + delta * real(npts-1)

  az = 0.0
  baz = 0.0
  dist = 0.0
  call distaz(evla,evlo,stla,stlo,dist,az,baz,gcarc,ier)

  if ( icomp==1 ) then
  cmpaz = 0.0
  cmpinc = 0.0
  else
  if ( icomp==2 ) then
  cmpaz = baz + 180.0
  if ( cmpaz>=360.0 ) cmpaz = cmpaz - 360.0
  cmpinc = 90.0
  else
  if ( icomp==3 ) then
  cmpaz = baz + 90.0
  if ( cmpaz>=360.0 ) cmpaz = cmpaz - 360.0
  cmpinc = 90.0
  else
  stop 'Illegal component type (wsact0)'
  endif
  endif
  endif

  nvhdr = 6
  ninf = 0
  nhst = 0
  iftype = 1
  leven = 1
  lovrok = 1
  lcalda = 1

  ier = filenchk( 80,output )
  ier = sacbw( mfhdr,mnhdr,mihdr,mlhdr,mkhdr, fhdr,nhdr,ihdr,lhdr,khdr, npts,y,output )

  end

!-----------------------------------------------------------------------

  subroutine fast_fourier_transform_real(input_values,nn,sign_of_transform)

  implicit none

  real(kind=8) input_values(*)
  integer nn,sign_of_transform
  real(kind=8) wr,wi,pr_ww,pi_ww,temporary_ww,angle
  real(kind=8) r_ww_temp,i_ww_temp
  integer n,i,j,m,maxmval,step_value

  real(kind=8), parameter :: two_pi = 6.28318530717959d0

  n=2*nn
  j=1

  do i=1,n,2
  if(j > i) then
    r_ww_temp=input_values(j)
    i_ww_temp=input_values(j+1)
    input_values(j)=input_values(i)
    input_values(j+1)=input_values(i+1)
    input_values(i)=r_ww_temp
    input_values(i+1)=i_ww_temp
  endif
  m=n/2
 1 if (m >= 2 .and. j > m) then
    j=j-m
    m=m/2
  goto 1
  endif
  j=j+m
  enddo

  maxmval=2
 2 if (n > maxmval) then
  step_value=2*maxmval
  angle = two_pi / (sign_of_transform*maxmval)
  pr_ww=-2.d0*dsin(0.5d0*angle)**2
  pi_ww=dsin(angle)
  wr=1.d0
  wi=0.d0

  do m=1,maxmval,2
    do i=m,n,step_value
      j=i+maxmval
      r_ww_temp=wr*input_values(j)-wi*input_values(j+1)
      i_ww_temp=wr*input_values(j+1)+wi*input_values(j)
      input_values(j)=input_values(i)-r_ww_temp
      input_values(j+1)=input_values(i+1)-i_ww_temp
      input_values(i)=input_values(i)+r_ww_temp
      input_values(i+1)=input_values(i+1)+i_ww_temp
    enddo
    temporary_ww=wr
    wr=wr*pr_ww-wi*pi_ww+wr
    wi=wi*pr_ww+temporary_ww*pi_ww+wi
  enddo

  maxmval=step_value
  goto 2
  endif

  end

