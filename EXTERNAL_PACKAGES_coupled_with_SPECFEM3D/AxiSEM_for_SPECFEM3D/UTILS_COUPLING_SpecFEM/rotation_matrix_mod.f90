  module rotation_matrix
  contains
    subroutine def_rot_matrix(srccolat,srclon,rot_mat,trans_rot_mat)

      use global_parameters, only: CUSTOM_REAL
      real(kind=CUSTOM_REAL), dimension(3,3)  :: rot_mat,trans_rot_mat
      integer                       :: i, j
      real(kind=CUSTOM_REAL)                  :: smallval,srccolat,srclon

      smallval=1e-11

      ! This is the rotation matrix of Nissen-Meyer, Dahlen, Fournier, GJI 2007.
      rot_mat(1,1) = dcos(srccolat) * dcos(srclon)
      rot_mat(2,2) = dcos(srclon)
      rot_mat(3,3) = dcos(srccolat)
      rot_mat(2,1) = dcos(srccolat) * dsin(srclon)
      rot_mat(3,1) = -dsin(srccolat)
      rot_mat(3,2) = 0.d0
      rot_mat(1,2) = -dsin(srclon)
      rot_mat(1,3) = dsin(srccolat) * dcos(srclon)
      rot_mat(2,3) = dsin(srccolat) * dsin(srclon)

      where (dabs(rot_mat) < smallval) rot_mat = 0.0

      trans_rot_mat = transpose(rot_mat)

      write(*,*) ' ROTATION MATRIX '
      write(*,*) rot_mat(1,1),rot_mat(1,2),rot_mat(1,3)
      write(*,*) rot_mat(2,1),rot_mat(2,2),rot_mat(2,3)
      write(*,*) rot_mat(3,1),rot_mat(3,2),rot_mat(3,3)
      write(*,*)
    end subroutine def_rot_matrix
! --------------------------------------------

subroutine def_rot_matrix_DG(srccolat,srclon,rot_mat,trans_rot_mat)

    use global_parameters, only: CUSTOM_REAL
    !integer, parameter :: CUSTOM_REAL=8
    real(kind=CUSTOM_REAL), dimension(3,3)  :: rot_mat,trans_rot_mat
    integer                       :: i, j
    real(kind=CUSTOM_REAL)                  :: smallval,srccolat,srclon

    smallval=1e-11

    ! This is the rotation matrix of Nissen-Meyer, Dahlen, Fournier, GJI 2007.
!!$    rot_mat(1,1) = -sin(srclon)
!!$    rot_mat(2,2) = cos(srccolat) * sin(srclon)
!!$    rot_mat(3,3) = cos(srccolat)
!!$
!!$    rot_mat(2,3) = sin(srccolat) * sin(srclon)
!!$    rot_mat(1,3) = sin(srccolat) * cos(srclon)
!!$    rot_mat(1,2) = cos(srccolat) * cos(srclon)
!!$
!!$    rot_mat(2,1) = cos(srclon)
!!$    rot_mat(3,2) = sin(srccolat)
!!$    rot_mat(3,1) = 0. !sin(srccolat) * sin(srclon)

!!$    rot_mat(1,1) = -cos(srccolat)*cos(srclon)
!!$    rot_mat(2,2) =  cos(srclon)
!!$    rot_mat(3,3) =  cos(srccolat)
!!$
!!$    rot_mat(2,3) = sin(srccolat) * sin(srclon)
!!$    rot_mat(1,3) = sin(srccolat) * cos(srclon)
!!$    rot_mat(1,2) = -sin(srclon)
!!$
!!$    rot_mat(2,1) = -cos(srccolat)*sin(srclon)
!!$    rot_mat(3,2) = 0.
!!$    rot_mat(3,1) = sin(srccolat) !sin(srccolat) * sin(srclon)

    rot_mat(1,2) = -cos(srccolat)*cos(srclon)
    rot_mat(2,1) =  cos(srclon)
    rot_mat(3,3) =  cos(srccolat)

    rot_mat(2,3) = sin(srccolat) * sin(srclon)
    rot_mat(1,3) = sin(srccolat) * cos(srclon)
    rot_mat(1,1) = -sin(srclon)

    rot_mat(2,2) = -cos(srccolat)*sin(srclon)
    rot_mat(3,1) = 0.
    rot_mat(3,2) = sin(srccolat) !sin(srccolat) * sin(srclon)


   where (abs(rot_mat) < smallval) rot_mat = 0.0

    write(*,*) rot_mat(1,1),rot_mat(1,2),rot_mat(1,3)
    write(*,*) rot_mat(2,1),rot_mat(2,2),rot_mat(2,3)
    write(*,*) rot_mat(3,1),rot_mat(3,2),rot_mat(3,3)

    trans_rot_mat = transpose(rot_mat)
  end subroutine def_rot_matrix_DG


!--------------------------------------------------

subroutine def_rot_azi_chunk(rot,trot,azi)

  use global_parameters, only: CUSTOM_REAL
  real(kind=CUSTOM_REAL)                  :: azi
  real(kind=CUSTOM_REAL), dimension(3,3)  :: rot,trot

  rot(1,1)=cos(azi)
  rot(1,2)=-1.0*sin(azi) !sin(azi)
  rot(1,3)=0._CUSTOM_REAL

  rot(2,1)=sin(azi) !-sin(azi)
  rot(2,2)=cos(azi)
  rot(2,3)=0._CUSTOM_REAL

  rot(3,1)=0._CUSTOM_REAL
  rot(3,2)=0._CUSTOM_REAL
  rot(3,3)=1._CUSTOM_REAL

  trot=transpose(rot)

end subroutine def_rot_azi_chunk

!----------------------------------------------------

  subroutine rotate_box(r,th,ph,trans_rot_mat)
    use global_parameters, only: CUSTOM_REAL

    real(kind=CUSTOM_REAL)    trans_rot_mat(3,3)
    real(kind=CUSTOM_REAL)   ,intent(inout) :: r,th,ph
    real(kind=CUSTOM_REAL)    :: x_vec(3), x_vec_rot(3), r_r,smallval_dble


    smallval_dble=0.d0 !! 1e-11  ! a quoi sert ce truc? ! VM VM
    !! verifier l'effet que ca peut avoir de le mettre a zero

    x_vec(1) = r * dsin(th) * dcos(ph)
    x_vec(2) = r * dsin(th) * dsin(ph)
    x_vec(3) = r * dcos(th)

    x_vec_rot = matmul(trans_rot_mat,x_vec)

    !write(23,*) x_vec
    !write(23,*) x_vec_rot

    r_r = dsqrt(x_vec_rot(1)**2 + x_vec_rot(2)**2 + x_vec_rot(3)**2 )
    th = dacos((x_vec_rot(3)  + smallval_dble )/ ( r_r + smallval_dble) )
    ph = atan2(x_vec_rot(2),x_vec_rot(1))

  end subroutine rotate_box

!----------------------------------------------
  end module rotation_matrix
