

     write(outputname,"('/scratch/DATABASES_MPI_DIMITRI/shaking_data',i8.8,'.sac')") ipoin

! create file, or append to end of file
     if(icounter_shakemap == 1) then
       open(unit=IOUT,file=outputname,status='unknown')
!! DK DK UGLY problem,  cannot know min max mean value in advance, therefore fictitious
       call write_sac2000_header(NSTEP,DT,-1.d0,+1.d0,0.d0,IOUT)
     else
       open(unit=IOUT,file=outputname,status='old',position='append')
     endif

! ZZZZZZZZZZZZZZZZZZZ   DK DK UGLY

! write data values by groups of NVALUES_SAC2000_LINE per line

  nb_samples = NSIZESAVE_SAC2000

! check if we need to write fewer samples for last time slice
  if((icounter_shakemap-1)*NSIZESAVE_SAC2000 + nb_samples > NSTEP) &
    nb_samples = NSTEP - (icounter_shakemap-1)*NSIZESAVE_SAC2000

  if(nb_samples <= 0) call exit_MPI(myrank,'error number of samples for last time slice')

  it_slice = 0

  do iline = 1,nb_samples/NVALUES_SAC2000_LINE
    write(IOUT,"(5g15.7)") (store_val_norm_veloc_timeslice(it_slice + ivalue,ipoin),ivalue=1,NVALUES_SAC2000_LINE)
    it_slice = it_slice + NVALUES_SAC2000_LINE
  enddo

! write remaining data values
  iremain = nb_samples - (nb_samples/NVALUES_SAC2000_LINE)*NVALUES_SAC2000_LINE
  if(iremain > 0) write(IOUT,"(5g15.7)") (store_val_norm_veloc_timeslice(it_slice + ivalue,ipoin),ivalue=1,iremain)

     close(IOUT)

    enddo   ! end of triple loop on points in surface elements
    enddo
    enddo


  end program specfem3D



!!!!!!!!!!!!!!!! DK DK UGLY     sac2000


  subroutine write_sac2000_header(NPTS,DELTA,DEPMIN,DEPMAX,DEPMEN,ifile)

! sac2000 alphanumeric format
! The header section is stored on the first 30 lines
! followed by the data section in 5g15.7 format

!! DK DK UGLY problem,  cannot know min max mean value in advance
!! DK DK UGLY but we can hope that they are not used (information only)

  implicit none

  include "constants.h"

! output file
  integer ifile

! number of time steps
  integer NPTS

! duration of time step
  double precision DELTA

! min, max and mean value of seismogram
  double precision DEPMIN,DEPMAX,DEPMEN

! a few particular constant values

! header version number
  integer, parameter :: NVHDR = 6

! flag to indicate time series file
  integer, parameter :: IFTYPE = 1

! flag to indicate that data is evenly spaced
  integer, parameter :: LEVEN = 1

! flag for station polarity
  integer, parameter :: LPSPOL = 0

! falg to indicate that it is okay to overwrite this file (not write-protected)
  integer, parameter :: LOVROK = 1

! flag for how to compute station coordinates
  integer, parameter :: LCALDA = 1

! unused values
  integer, parameter :: unused_int_null = 0
  integer, parameter :: unused_int_regular = -12345
  double precision, parameter :: internal_dble = -12345.d0
  double precision, parameter :: unused_dble = -12345.d0

  integer NZYEAR,NZJDAY,NZHOUR,NZMIN,NZSEC,NZMSEC,NORID,NEVID
  integer NSPTS,NWFID,NXSIZE,NYSIZE,IDEP,IZTYPE,IINST
  integer ISTREG,IEVREG,IEVTYP,IQUAL,ISYNTH,IMAGTYP,IMAGSRC

  double precision SCALE,ODELTA,B,E,O,A
  double precision T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,F,RESP0,RESP1,RESP2,RESP3
  double precision RESP4,RESP5,RESP6,RESP7,RESP8,RESP9,STLA,STLO,STEL,STDP
  double precision EVLA,EVLO,EVEL,EVDP,MAG,USER0,USER1,USER2,USER3,USER4
  double precision USER5,USER6,USER7,USER8,USER9,DIST,AZ,BAZ,GCARC
  double precision CMPAZ,CMPINC,XMINIMUM
  double precision XMAXIMUM,YMINIMUM,YMAXIMUM,ADJTM

! absolute time of first and last time steps
!! DK DK UGLY this is slightly wrong strictly speaking, use -hdur /2
  B = 0.d0
  E = (NPTS-1) * DELTA

! define default values for other variables (always -12345 in sac2000)
  NZYEAR = -12345
  NZJDAY = -12345
  NZHOUR = -12345
  NZMIN = -12345
  NZSEC = -12345
  NZMSEC = -12345
  NORID = -12345
  NEVID = -12345
  NSPTS = -12345
  NWFID = -12345
  NXSIZE = -12345
  NYSIZE = -12345
  IDEP = -12345
  IZTYPE = -12345
  IINST = -12345
  ISTREG = -12345
  IEVREG = -12345
  IEVTYP = -12345
  IQUAL = -12345
  ISYNTH = -12345
  IMAGTYP = -12345
  IMAGSRC = -12345

  SCALE = -12345.d0
  ODELTA = -12345.d0
  O = -12345.d0
  A = -12345.d0
  T0 = -12345.d0
  T1 = -12345.d0
  T2 = -12345.d0
  T3 = -12345.d0
  T4 = -12345.d0
  T5 = -12345.d0
  T6 = -12345.d0
  T7 = -12345.d0
  T8 = -12345.d0
  T9 = -12345.d0
  F = -12345.d0
  RESP0 = -12345.d0
  RESP1 = -12345.d0
  RESP2 = -12345.d0
  RESP3 = -12345.d0
  RESP4 = -12345.d0
  RESP5 = -12345.d0
  RESP6 = -12345.d0
  RESP7 = -12345.d0
  RESP8 = -12345.d0
  RESP9 = -12345.d0
  STLA = -12345.d0
  STLO = -12345.d0
  STEL = -12345.d0
  STDP = -12345.d0
  EVLA = -12345.d0
  EVLO = -12345.d0
  EVEL = -12345.d0
  EVDP = -12345.d0
  MAG = -12345.d0
  USER0 = -12345.d0
  USER1 = -12345.d0
  USER2 = -12345.d0
  USER3 = -12345.d0
  USER4 = -12345.d0
  USER5 = -12345.d0
  USER6 = -12345.d0
  USER7 = -12345.d0
  USER8 = -12345.d0
  USER9 = -12345.d0
  DIST = -12345.d0
  AZ = -12345.d0
  BAZ = -12345.d0
  GCARC = -12345.d0
  CMPAZ = -12345.d0
  CMPINC = -12345.d0
  XMINIMUM = -12345.d0
  XMAXIMUM = -12345.d0
  YMINIMUM = -12345.d0
  YMAXIMUM = -12345.d0
  ADJTM = -12345.d0

! header section with reals

! header line 01
  write(ifile,"(5g15.7)") DELTA,DEPMIN,DEPMAX,SCALE,ODELTA

! header line 02
  write(ifile,"(5g15.7)") B,E,O,A,internal_dble

! header line 03
  write(ifile,"(5g15.7)") T0,T1,T2,T3,T4

! header line 04
  write(ifile,"(5g15.7)") T5,T6,T7,T8,T9

! header line 05
  write(ifile,"(5g15.7)") F,RESP0,RESP1,RESP2,RESP3

! header line 06
  write(ifile,"(5g15.7)") RESP4,RESP5,RESP6,RESP7,RESP8

! header line 07
  write(ifile,"(5g15.7)") RESP9,STLA,STLO,STEL,STDP

! header line 08
  write(ifile,"(5g15.7)") EVLA,EVLO,EVEL,EVDP,MAG

! header line 09
  write(ifile,"(5g15.7)") USER0,USER1,USER2,USER3,USER4

! header line 10
  write(ifile,"(5g15.7)") USER5,USER6,USER7,USER8,USER9

! header line 11
  write(ifile,"(5g15.7)") DIST,AZ,BAZ,GCARC,internal_dble

! header line 12
  write(ifile,"(5g15.7)") internal_dble,DEPMEN,CMPAZ,CMPINC,XMINIMUM

! header line 13
  write(ifile,"(5g15.7)") XMAXIMUM,YMINIMUM,YMAXIMUM,ADJTM,unused_dble

! header line 14
  write(ifile,"(5g15.7)") unused_dble,unused_dble,unused_dble,unused_dble,unused_dble


! header section with integers

! header line 15
  write(ifile,"(5i10)") NZYEAR,NZJDAY,NZHOUR,NZMIN,NZSEC

! header line 16
  write(ifile,"(5i10)") NZMSEC,NVHDR,NORID,NEVID,NPTS

! header line 17
  write(ifile,"(5i10)") NSPTS,NWFID,NXSIZE,NYSIZE,unused_int_regular

! header line 18
  write(ifile,"(5i10)") IFTYPE,IDEP,IZTYPE,unused_int_regular,IINST

! header line 19
  write(ifile,"(5i10)") ISTREG,IEVREG,IEVTYP,IQUAL,ISYNTH

! header line 20
  write(ifile,"(5i10)") IMAGTYP,IMAGSRC,unused_int_regular,unused_int_regular,unused_int_regular

! header line 21
  write(ifile,"(5i10)") unused_int_regular,unused_int_regular,unused_int_regular,unused_int_regular,unused_int_regular

! header line 22
  write(ifile,"(5i10)") LEVEN,LPSPOL,LOVROK,LCALDA


! header section with strings

! header line 23
  write(ifile,"(a25)") '-12345           -12345  '

! header line 24
  write(ifile,"(a24)") '-12345  -12345  -12345  '

! header line 25
  write(ifile,"(a24)") '-12345  -12345  -12345  '

! header line 26
  write(ifile,"(a24)") '-12345  -12345  -12345  '

! header line 27
  write(ifile,"(a24)") '-12345  -12345  -12345  '

! header line 28
  write(ifile,"(a24)") '-12345  -12345  -12345  '

! header line 29
  write(ifile,"(a24)") '-12345  -12345  -12345  '

! header line 30
  write(ifile,"(a24)") '-12345  -12345  -12345  '

  end subroutine write_sac2000_header

