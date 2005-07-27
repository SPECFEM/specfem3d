!
! Dimitri Komatitsch, University of Pau, France, April 2004
!
!  rotate c_ijkl for copper crystal study with Jose Carcione
!
! program to rotate the 21 components of the full 3D anisotropic c_ijkl tensor
!
! based on Helbig, Foundations of anisotropy for exploration seismics, 1994
!
! using the formulas of p. 72 and the listing of appendix 3C p. 123
!
  program rotate_c_ijkl_helbig

  implicit none

! axis and angle (in degree) of rotation
  integer, parameter :: iaxis = 3
  double precision, parameter :: angle = 45.d0

!! DK DK
!! DK DK copper crystal (cubic symmetry)
!! DK DK
  double precision, parameter :: c11 = 169.d9
  double precision, parameter :: c12 = 122.d9
  double precision, parameter :: c13 = c12
  double precision, parameter :: c14 = 0.d0
  double precision, parameter :: c15 = 0.d0
  double precision, parameter :: c16 = 0.d0

  double precision, parameter :: c22 = c11
  double precision, parameter :: c23 = c12
  double precision, parameter :: c24 = 0.d0
  double precision, parameter :: c25 = 0.d0
  double precision, parameter :: c26 = 0.d0

  double precision, parameter :: c33 = c11
  double precision, parameter :: c34 = 0.d0
  double precision, parameter :: c35 = 0.d0
  double precision, parameter :: c36 = 0.d0

  double precision, parameter :: c44 = 75.3d9
  double precision, parameter :: c45 = 0.d0
  double precision, parameter :: c46 = 0.d0

  double precision, parameter :: c55 = c44
  double precision, parameter :: c56 = 0.d0

  double precision, parameter :: c66 = c44

! reduced c_IJ parameters in input, before rotation
  double precision c_in(6,6)

! reduced c_IJ parameters in output, after rotation
  double precision c_out(6,6)

! full c_ijkl parameters in input, before rotation
  double precision c_full_in(3,3,3,3)

! full c_ijkl parameters in output, after rotation
  double precision c_full_out(3,3,3,3)

  integer i,j

! assign initial value to the input matrix in reduced notation
  c_in(1,1) = c11
  c_in(2,2) = c22
  c_in(3,3) = c33
  c_in(4,4) = c44
  c_in(5,5) = c55
  c_in(6,6) = c66
  c_in(1,2) = c12
  c_in(1,3) = c13
  c_in(2,3) = c23
  c_in(1,4) = c14
  c_in(2,4) = c24
  c_in(3,4) = c34
  c_in(1,5) = c15
  c_in(2,5) = c25
  c_in(3,5) = c35
  c_in(4,5) = c45
  c_in(1,6) = c16
  c_in(2,6) = c26
  c_in(3,6) = c36
  c_in(4,6) = c46
  c_in(5,6) = c56

! implement the symmetry of the reduced matrix
  do i=1,6
    do j=1,6
      if(j < i) c_in(i,j) = c_in(j,i)
    enddo
  enddo

! convert stiffnesses c_IJ to c_ijkl (reduced to full notation)
! using Helbig's routine
  call stiffness(c_full_in,c_in,1)

! rotate the full c_ijkl tensor
  call rotate4(c_full_in,c_full_out,iaxis,angle)

! convert stiffnesses c_ijkl to c_IJ (full to reduced notation)
! using Helbig's routine
  call stiffness(c_full_out,c_out,2)

! output the symmetric part of the reduced matrix before rotation
  print *
  print *,'! original c_ijkl before rotation'
  print *
  do i=1,6
    do j=1,6
      if(j >= i) write(*,100) i,j,c_in(i,j)
    enddo
  enddo
  print *

! output the symmetric part of the reduced rotated matrix
  print *
  print *,'! new c_ijkl after rotation'
  print *
  do i=1,6
    do j=1,6
      if(j >= i) write(*,100) i,j,c_out(i,j)
    enddo
  enddo
  print *

 100 format('  double precision, parameter :: c',i1,i1,'_copper = ',e30.20)

  end program rotate_c_ijkl_helbig


!----
!---- Helbig's code to convert c_IJ to c_ijkl (reduced to full notation)
!---- if mode=1, or c_ijkl to c_IJ (full to reduced notation) if mode=2
!----

  subroutine stiffness(a4,b2,mode)

  implicit none

  double precision a4(3,3,3,3),b2(6,6)
  integer mode

  integer i,j,k,l,p,q

  do i=1,3
    do j=1,3

      if(i == j) then
        p=i
      else
        p=9-i-j
      endif

      do k=1,3
        do l=1,3
          if(k == l) then
            q=k
          else
            q=9-k-l
          endif

          if(mode == 1) then
            a4(i,j,k,l) = b2(p,q)
          else if(mode == 2) then
            b2(p,q) = a4(i,j,k,l)
          else
            stop 'improper mode'
          endif

        enddo
      enddo

    enddo
  enddo

  end subroutine stiffness

!----
!---- Helbig's code to rotate the full c_ijkl
!----

  subroutine rotate4(c,cs,iax,alfa)

! c,cs: pre-rotation and post-rotation stiffnesses
! iax: axis of rotation (between 1 and 3)
! alfa: counterclockwise angle in degrees
! the routine moves the 1-axis towards the 2-axis if iax=3

  implicit none

  double precision c(3,3,3,3),cs(3,3,3,3),a(3,3)
  double precision alfa,al,cc,ss
  integer i,j,k,l,is,js,ks,ls,iax,i1,i2,i3

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0

  if(iax < 1 .or. iax > 3) stop 'wrong axis number'

! select axes
  i3=iax
  i1=i3+1
  if(i1 > 3) i1=i1-3
  i2=i1+1
  if(i2 > 3) i2=i2-3

! set up transformation matrix
  al = alfa * PI/180.d0
  cc = dcos(al)
  ss = dsin(al)
  a(i1,i1) = cc
  a(i1,i2) = - ss
  a(i1,i3) = 0.d0
  a(i2,i1) = ss
  a(i2,i2) = cc
  a(i2,i3) = 0.d0
  a(i3,i1) = 0.d0
  a(i3,i2) = 0.d0
  a(i3,i3) = 1.d0

! actual rotation using temporary arrays
  do is=1,3
    do js=1,3
      do ks=1,3
        do ls=1,3
          cs(is,js,ks,ls) = 0.d0
          do i=1,3
            do j=1,3
              do k=1,3
                do l=1,3
      cs(is,js,ks,ls) = cs(is,js,ks,ls) + c(i,j,k,l)*a(i,is)*a(j,js)*a(k,ks)*a(l,ls)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  end subroutine rotate4

