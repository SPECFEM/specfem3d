module cmt3d_sub4

  use cmt3d_constants
  implicit none

contains

  !==================================================================

  subroutine rotate_cmt(cmt_par,npar,elat0,elon0,DIRECTION)

    ! DIRECTION = 1/-1: entering/leaving with cmt par read from /write to cmtsolution file
    ! in terms of Mij, dep, lon and lat

    real*8, intent(inout) :: cmt_par(NPARMAX)
    real*8, intent(in) :: elat0, elon0
    integer,intent(in)  :: npar,DIRECTION

    real*8 :: th, phi, loc(3), gl(3), moment(3,3), elon,elat,edep
    real*8 :: rmat(3,3)

    ! check input arguments
    if (DIRECTION /= 1 .and. DIRECTION /= -1) stop 'Error DIRECTION (1 or -1)'

    if (DIRECTION == 1) then ! local to global coordinates (dep,lon/E,lat/N) - > (r, th, phi)

       ! rmat converts local {r,t,p} (or {Z,S,E}) coord. to global {X,Y,Z} coord.
       ! X_global= rmat * X_local, then M_global= R * M_local * R^T
       call calc_rot_matrix(elon0,elat0,rmat)

       ! moment tensors
       moment(1,1)=cmt_par(1); moment(2,2)=cmt_par(2); moment(3,3)=cmt_par(3)
       moment(1,2)=cmt_par(4); moment(1,3)=cmt_par(5); moment(2,3)=cmt_par(6)
       moment(2,1)=moment(1,2); moment(3,1)=moment(1,3); moment(3,2)=moment(2,3)
       moment=matmul(rmat, matmul(moment,transpose(rmat)))
       cmt_par(1)=moment(1,1); cmt_par(2)=moment(2,2); cmt_par(3)=moment(3,3)
       cmt_par(4)=moment(1,2); cmt_par(5)=moment(1,3); cmt_par(6)=moment(2,3)

       if (npar > 6) then ! guaranteed to be npar >= 9
          edep=cmt_par(7); elon=cmt_par(8); elat=cmt_par(9)
          if (abs(elon-elon0) > EPS5 .or. abs(elat-elat0) > EPS5) stop 'Error lat,lon'
          loc=(/R_EARTH-edep,0.d1,0.d1/)  ! local [R,T,P]
          loc=matmul(rmat,loc)  ! global [X,Y,Z]
          cmt_par(7:9)=loc(1:3)  ![R,T,P] in local or [X,Y,Z] in global
       endif


    else ! DIRECTION = -1 : global to local coordinates

       if (npar > 6) then
          gl=cmt_par(7:9)/sqrt(sum(cmt_par(7:9)**2)) ! global [X,Y,Z]
          th=acos(gl(3))
          if (abs(th) < EPS5 .or. abs(th-pi) < EPS5) then
             phi=0
          else
             phi=atan2(gl(2)/sin(th),gl(1)/sin(th))
          endif
          elon=phi*180/PI; elat=90-th*180/PI
          if (elon > 180) elon=elon-360
       else
          elat=elat0; elon=elon0
       endif
       call calc_rot_matrix(elon,elat,rmat)

       ! moment tensor elements
       moment(1,1)=cmt_par(1); moment(2,2)=cmt_par(2); moment(3,3)=cmt_par(3)
       moment(1,2)=cmt_par(4); moment(1,3)=cmt_par(5); moment(2,3)=cmt_par(6)
       moment(2,1)=moment(1,2); moment(3,1)=moment(1,3); moment(3,2)=moment(2,3)
       moment=matmul(transpose(rmat), matmul(moment,rmat))
       cmt_par(1)=moment(1,1); cmt_par(2)=moment(2,2); cmt_par(3)=moment(3,3)
       cmt_par(4)=moment(1,2); cmt_par(5)=moment(1,3); cmt_par(6)=moment(2,3)

       if (npar > 6) then
          loc=matmul(transpose(rmat),cmt_par(7:9))  ! moment, loc in [R,T,P] local coord
          if (abs(loc(2)) > EPS5 .or. abs(loc(3)) > EPS5) stop 'Error loc(2,3)'
          cmt_par(7)=R_EARTH-loc(1); cmt_par(8)=elon; cmt_par(9)=elat  ! [depth,lon,lat]
       endif
    endif

  end subroutine rotate_cmt

  !==================================================================

  subroutine calc_rot_matrix(elon,elat,rmat)

    real*8,intent(in) :: elon, elat
    real*8,intent(out) :: rmat(3,3)
    real*8 :: th,phi

    th=(90-elat)*PI/180
    phi=elon*PI/180

    rmat(1,1)=dsin(th)*dcos(phi)
    rmat(1,2)=dcos(th)*dcos(phi)
    rmat(1,3)=-dsin(phi)
    rmat(2,1)=dsin(th)*dsin(phi)
    rmat(2,2)=dcos(th)*dsin(phi)
    rmat(2,3)=dcos(phi)
    rmat(3,1)=dcos(th)
    rmat(3,2)=-dsin(th)
    rmat(3,3)=0

  end subroutine calc_rot_matrix

end module cmt3d_sub4
