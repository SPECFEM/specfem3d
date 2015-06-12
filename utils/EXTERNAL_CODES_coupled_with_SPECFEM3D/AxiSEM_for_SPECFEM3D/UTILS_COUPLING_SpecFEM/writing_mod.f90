module writing_mod

  contains
  subroutine write_veloc3D(ivx,ivy,ivz)
    use global_parameters,only : nbrec,data_rec
    integer ivx,ivy,ivz
!!$    write(ivx) data_rec(:,1)
!!$    write(ivy) data_rec(:,2)
!!$    write(ivz) data_rec(:,3)

    ! permutation 
    write(ivx) data_rec(:,1)
    write(ivy) data_rec(:,2)
    write(ivz) data_rec(:,3)

    !write(*,*) ' write veloc ',data_rec(10,3)

  end subroutine write_veloc3D

  subroutine write_stress3D(isxx,isyy,iszz,isxy,isxz,isyz)
    use global_parameters,only :nbrec,stress_rec,SINGLE_REAL,mat,tmat,strain_rec
    integer isxx,isyy,iszz,isxy,isxz,isyz
    integer i,j,k,irec
    real(kind=SINGLE_REAL) tmp(3,3),tmp1(3,3),st(3,3)
    real(kind=SINGLE_REAL) lam,mu

!!$    do irec=1,nbrec
!!$           
!!$       ! stress cartesian
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
!!$             end do
!!$          end do
!!$       end do
!!$       
!!$       ! R*st*Rt
!!$       tmp1=0.
!!$       do j=1,3
!!$          do i=1,3
!!$             do k=1,3
!!$                tmp1(i,j)=tmp1(i,j)+tmp(i,k)*tmat(k,j)
!!$             end do
!!$          end do
!!$       end do
!!$       
!!$       ! stress in cartesian
!!$       stress_rec(irec,1)=tmp1(1,1)
!!$       stress_rec(irec,2)=tmp1(2,2)
!!$       stress_rec(irec,3)=tmp1(3,3)
!!$       stress_rec(irec,4)=tmp1(1,2) 
!!$       stress_rec(irec,5)=tmp1(1,3)
!!$       stress_rec(irec,6)=tmp1(2,3)
!!$       
!!$    end do

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

    !write(*,*)  stress_rec(:,1)
    
  end subroutine write_stress3D
end module writing_mod
