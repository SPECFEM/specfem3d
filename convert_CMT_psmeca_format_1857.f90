
! DK DK convert CMTSOLUTION to psmeca format to check finite sources for L.A.

  program convert

  implicit none

! number of lines in CMTSOLUTION file
!  integer, parameter :: NLINES = 14118
  integer, parameter :: NLINES = 3000

  integer iline
  double precision long,lat,depth,mrr,mtt,mpp,mrt,mrp,mtp

  character(len=150) string

! header of script
  write(*,8)

  write(*,*) 'psbasemap -R-124/-114/30/40 -JM15c -Bf1a2:Distance:/:"samples":WeSn -K > f1.ps'

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

! use scale of 1.e14

 Mrr = Mrr / 1.d14
 Mtt = Mtt / 1.d14
 Mpp = Mpp / 1.d14
 Mrt = Mrt / 1.d14
 Mrp = Mrp / 1.d14
 Mtp = Mtp / 1.d14

!! DK DK need to remove the last -K in output file
 if(iline < NLINES) then
   write(*,10) long,lat,depth,mrr,mtt,mpp,mrt,mrp,mtp
 else
   write(*,11) long,lat,depth,mrr,mtt,mpp,mrt,mrp,mtp
 endif

 8 format ('#!/bin/csh')

 10 format ('echo "',f10.4,' ',f10.4,' ',f10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4, &
   ' 14 0 0 " | psmeca -R-124/-114/30/40 -JM15c -Sm0.7c -O -K -G200/0/0 >> f1.ps')

 11 format ('echo "',f10.4,' ',f10.4,' ',f10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4,' ',e10.4, &
   ' 14 0 0 " | psmeca -R-124/-114/30/40 -JM15c -Sm0.7c -O -G200/0/0 >> f1.ps')

  enddo

  end

