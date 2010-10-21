
! DK DK convert CMTSOLUTION to psmeca format to check finite sources for L.A.

  program convert_CMT_psmeca_format

  implicit none

! number of lines in CMTSOLUTION file

!! DK DK choose between Northridge and SAF 1857
  integer, parameter :: SAF_1857 = 1
  integer, parameter :: NORTHRIDGE = 2
  integer, parameter :: ichoice = NORTHRIDGE

  integer iline,NLINES
  double precision long,lat,depth,mrr,mtt,mpp,mrt,mrp,mtp,scaleval

  character(len=256) string

! header of script
  write(*,8)

  if(ichoice == SAF_1857) then
!   NLINES = 14118   !! DK DK real value
    NLINES = 3000    !! DK DK use only part of the file
    scaleval = 1.d14
    write(*,*) 'psbasemap -R-124/-114/30/40 -JM15c -Bf1a2:Distance:/:"samples":WeSn -K > f1.ps'
  else if(ichoice == NORTHRIDGE) then
    NLINES = 196
    scaleval = 1.d24
    write(*,*) 'psbasemap -R-120/-116/32/36 -JM15c -Bf1a2:Distance:/:"samples":WeSn -K > f1.ps'
  else
    stop 'incorrect value of ichoice'
  endif

! Format for psmeca: lon lat depth mrr mtt mpp mrt mrp mtp iexp name

  do iline = 1,NLINES

! skip four lines of labels
    read(*,*)
    read(*,*)
    read(*,*)
    read(*,*)

! read latitude
  read(*,"(a)") string
  read(string(10:len_trim(string)),*) lat

! read longitude
  read(*,"(a)") string
  read(string(11:len_trim(string)),*) long

! read depth
  read(*,"(a)") string
  read(string(7:len_trim(string)),*) depth

! read Mrr
  read(*,"(a)") string
  read(string(5:len_trim(string)),*) Mrr

! read Mtt
  read(*,"(a)") string
  read(string(5:len_trim(string)),*) Mtt

! read Mpp
  read(*,"(a)") string
  read(string(5:len_trim(string)),*) Mpp

! read Mrt
  read(*,"(a)") string
  read(string(5:len_trim(string)),*) Mrt

! read Mrp
  read(*,"(a)") string
  read(string(5:len_trim(string)),*) Mrp

! read Mtp
  read(*,"(a)") string
  read(string(5:len_trim(string)),*) Mtp

! use normalized scale
  Mrr = Mrr / scaleval
  Mtt = Mtt / scaleval
  Mpp = Mpp / scaleval
  Mrt = Mrt / scaleval
  Mrp = Mrp / scaleval
  Mtp = Mtp / scaleval

!! DK DK need to remove the last -K in output file
  if(iline < NLINES) then
    if(ichoice == SAF_1857) then
      write(*,100) long,lat,depth,mrr,mtt,mpp,mrt,mrp,mtp
    else
      write(*,200) long,lat,depth,mrr,mtt,mpp,mrt,mrp,mtp
    endif
  else
    if(ichoice == SAF_1857) then
      write(*,110) long,lat,depth,mrr,mtt,mpp,mrt,mrp,mtp
    else
      write(*,210) long,lat,depth,mrr,mtt,mpp,mrt,mrp,mtp
    endif
  endif

  enddo

 8 format ('#!/bin/csh')

 100 format ('echo "',f10.4,' ',f10.4,' ',f10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4, &
   ' 14 0 0 " | psmeca -R-124/-114/30/40 -JM15c -Sm0.7c -O -K -G200/0/0 >> f1.ps')

 110 format ('echo "',f10.4,' ',f10.4,' ',f10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4, &
   ' 14 0 0 " | psmeca -R-124/-114/30/40 -JM15c -Sm0.7c -O -G200/0/0 >> f1.ps')

 200 format ('echo "',f10.4,' ',f10.4,' ',f10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4, &
   ' 24 0 0 " | psmeca -R-120/-116/32/36 -JM15c -Sm0.7c -O -K -G200/0/0 >> f1.ps')

 210 format ('echo "',f10.4,' ',f10.4,' ',f10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4, &
   ' 24 0 0 " | psmeca -R-120/-116/32/36 -JM15c -Sm0.7c -O -G200/0/0 >> f1.ps')

  end program convert_CMT_psmeca_format

