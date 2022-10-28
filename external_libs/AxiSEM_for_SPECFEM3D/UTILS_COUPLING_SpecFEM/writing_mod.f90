module writing_mod

  contains

!
!=====================================================================================!
!

  subroutine write_deriv_3D(idu1d1,idu1d2,idu1d3,idu2d1,idu2d2,idu2d3,idu3d1,idu3d2,idu3d3)

    use global_parameters, only: nbrec,deriv_rec
    integer idu1d1,idu1d2,idu1d3,idu2d1,idu2d2,idu2d3,idu3d1,idu3d2,idu3d3

    write(idu1d1) deriv_rec(:,1)
    write(idu1d2) deriv_rec(:,2)
    write(idu1d3) deriv_rec(:,3)

    write(idu2d1) deriv_rec(:,4)
    write(idu2d2) deriv_rec(:,5)
    write(idu2d3) deriv_rec(:,6)

    write(idu3d1) deriv_rec(:,7)
    write(idu3d2) deriv_rec(:,8)
    write(idu3d3) deriv_rec(:,9)


!!$    write(idu1d1,'(2x, 7e20.10)') deriv_rec(1,1), deriv_rec(2,1), deriv_rec(3,1), &
!!$                                  deriv_rec(4,1), deriv_rec(5,1), deriv_rec(6,1), deriv_rec(7,1)
!!$
!!$    write(idu1d2,'(2x, 7e20.10)') deriv_rec(1,2), deriv_rec(2,2), deriv_rec(3,2), &
!!$                                  deriv_rec(4,2), deriv_rec(5,2), deriv_rec(6,2), deriv_rec(7,2)
!!$
!!$    write(idu1d3,'(2x, 7e20.10)') deriv_rec(1,3), deriv_rec(2,3), deriv_rec(3,3), &
!!$                                  deriv_rec(4,3), deriv_rec(5,3), deriv_rec(6,3), deriv_rec(7,3)
!!$
!!$    write(idu2d1,'(2x, 7e20.10)') deriv_rec(1,4), deriv_rec(2,4), deriv_rec(3,4), &
!!$                                  deriv_rec(4,4), deriv_rec(5,4), deriv_rec(6,4), deriv_rec(7,4)
!!$
!!$    write(idu2d2,'(2x, 7e20.10)') deriv_rec(1,5), deriv_rec(2,5), deriv_rec(3,5), &
!!$                                  deriv_rec(4,5), deriv_rec(5,5), deriv_rec(6,5), deriv_rec(7,5)
!!$
!!$    write(idu2d3,'(2x, 7e20.10)') deriv_rec(1,6), deriv_rec(2,6), deriv_rec(3,6), &
!!$                                  deriv_rec(4,6), deriv_rec(5,6), deriv_rec(6,6), deriv_rec(7,6)
!!$
!!$    write(idu3d1,'(2x, 7e20.10)') deriv_rec(1,7), deriv_rec(2,7), deriv_rec(3,7), &
!!$                                  deriv_rec(4,7), deriv_rec(5,7), deriv_rec(6,7), deriv_rec(7,7)
!!$
!!$    write(idu3d2,'(2x, 7e20.10)') deriv_rec(1,8), deriv_rec(2,8), deriv_rec(3,8), &
!!$                                  deriv_rec(4,8), deriv_rec(5,8), deriv_rec(6,8), deriv_rec(7,8)
!!$
!!$    write(idu3d3,'(2x, 7e20.10)') deriv_rec(1,9), deriv_rec(2,9), deriv_rec(3,9), &
!!$                                  deriv_rec(4,9), deriv_rec(5,9), deriv_rec(6,9), deriv_rec(7,9)


  end subroutine write_deriv_3D

!
!=====================================================================================!
!

  subroutine write_veloc_or_displ_3D(ivx,ivy,ivz)

    use global_parameters, only: nbrec,data_rec
    integer ivx,ivy,ivz
!!$    write(ivx) data_rec(:,1)
!!$    write(ivy) data_rec(:,2)
!!$    write(ivz) data_rec(:,3)


    write(ivx) data_rec(:,1)
    write(ivy) data_rec(:,2)
    write(ivz) data_rec(:,3)


!!$    write(ivx,'(2x, 7e20.10)') data_rec(1,1), data_rec(2,1), data_rec(3,1), &
!!$                               data_rec(4,1), data_rec(5,1), data_rec(6,1), data_rec(7,1)
!!$
!!$    write(ivy,'(2x, 7e20.10)') data_rec(1,2), data_rec(2,2), data_rec(3,2), &
!!$                               data_rec(4,2), data_rec(5,2), data_rec(6,2), data_rec(7,2)
!!$
!!$    write(ivz,'(2x, 7e20.10)') data_rec(1,3), data_rec(2,3), data_rec(3,3), &
!!$                               data_rec(4,3), data_rec(5,3), data_rec(6,3), data_rec(7,3)
!!$
!!$    !write(*,*) ' write veloc ',data_rec(10,3)


  end subroutine write_veloc_or_displ_3D

!
!=====================================================================================!
!

  subroutine write_stress3D(isxx,isyy,iszz,isxy,isxz,isyz)
    use global_parameters, only: nbrec,stress_rec,SINGLE_REAL,mat,tmat,strain_rec
    integer isxx,isyy,iszz,isxy,isxz,isyz
    integer i,j,k,irec
    real(kind=SINGLE_REAL) tmp(3,3),tmp1(3,3),st(3,3)
    real(kind=SINGLE_REAL) lam,mu

!!$    do irec=1,nbrec
!!$
!!$       ! stress Cartesian
!!$       st(1,1)=stress_rec(irec,1)
!!$       st(1,2)=stress_rec(irec,4)
!!$       st(1,3)=stress_rec(irec,5)
!!$
!!$       st(2,1)=stress_rec(irec,4)
!!$       st(2,2)=stress_rec(irec,2)
!!$       st(2,3)=stress_rec(irec,6)
!!$
!!$       st(3,1)=stress_rec(irec,5)
!!$       st(3,2)=stress_rec(irec,6)
!!$       st(3,3)=stress_rec(irec,3)
!!$
!!$       !
!!$       ! R*st
!!$       tmp=0.
!!$       do j=1,3
!!$          do i=1,3
!!$             do k=1,3
!!$                tmp(i,j)=tmp(i,j)+mat(i,k)*st(k,j)
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$       ! R*st*Rt
!!$       tmp1=0.
!!$       do j=1,3
!!$          do i=1,3
!!$             do k=1,3
!!$                tmp1(i,j)=tmp1(i,j)+tmp(i,k)*tmat(k,j)
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$       ! stress in Cartesian
!!$       stress_rec(irec,1)=tmp1(1,1)
!!$       stress_rec(irec,2)=tmp1(2,2)
!!$       stress_rec(irec,3)=tmp1(3,3)
!!$       stress_rec(irec,4)=tmp1(1,2)
!!$       stress_rec(irec,5)=tmp1(1,3)
!!$       stress_rec(irec,6)=tmp1(2,3)
!!$
!!$    enddo

    !strain_rec=stress_rec

    !lam=2.205e+11
    !mu=1.6200e+11

    !stress_rec(:,1)=(lam+2.*mu)*(strain_rec(:,1) + strain_rec(:,2))-2.*mu*strain_rec(:,2)+lam*strain_rec(:,3)
    !stress_rec(:,2)=(lam+2.*mu)*(strain_rec(:,1) + strain_rec(:,2))-2.*mu*strain_rec(:,1)+lam*strain_rec(:,3)
    !stress_rec(:,3)=lam*(strain_rec(:,1) + strain_rec(:,2)) + (lam+2.*mu)*strain_rec(:,3)
    !stress_rec(:,4)=2*mu*strain_rec(:,4)
    !stress_rec(:,5)=2*mu*strain_rec(:,5)
    !stress_rec(:,6)=2*mu*strain_rec(:,6)


    write(isxx) stress_rec(:,1)
    write(isyy) stress_rec(:,2)
    write(iszz) stress_rec(:,3)
    write(isxy) stress_rec(:,4)
    write(isxz) stress_rec(:,5)
    write(isyz) stress_rec(:,6)


!!$    write(isxx,'(2x, 7e20.10)') stress_rec(1,1), stress_rec(2,1), stress_rec(3,1), &
!!$                                  stress_rec(4,1), stress_rec(5,1), stress_rec(6,1), stress_rec(7,1)
!!$
!!$    write(isyy,'(2x, 7e20.10)') stress_rec(1,2), stress_rec(2,2), stress_rec(3,2), &
!!$                                  stress_rec(4,2), stress_rec(5,2), stress_rec(6,2), stress_rec(7,2)
!!$
!!$    write(iszz,'(2x, 7e20.10)') stress_rec(1,3), stress_rec(2,3), stress_rec(3,3), &
!!$                                  stress_rec(4,3), stress_rec(5,3), stress_rec(6,3), stress_rec(7,3)
!!$
!!$    write(isxy,'(2x, 7e20.10)') stress_rec(1,4), stress_rec(2,4), stress_rec(3,4), &
!!$                                  stress_rec(4,4), stress_rec(5,4), stress_rec(6,4), stress_rec(7,4)
!!$
!!$    write(isxz,'(2x, 7e20.10)') stress_rec(1,5), stress_rec(2,5), stress_rec(3,5), &
!!$                                  stress_rec(4,5), stress_rec(5,5), stress_rec(6,5), stress_rec(7,5)
!!$
!!$    write(isyz,'(2x, 7e20.10)') stress_rec(1,6), stress_rec(2,6), stress_rec(3,6), &
!!$                                  stress_rec(4,6), stress_rec(5,6), stress_rec(6,6), stress_rec(7,6)
!!$
!!$
!!$    !write(*,*)  stress_rec(:,1)


  end subroutine write_stress3D
end module writing_mod
