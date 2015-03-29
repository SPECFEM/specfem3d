module deriv_DF_mod
  include 'constants.h'
contains



  subroutine sub_D1(f,f1,nx,ny,nz,h)

    implicit none
    integer nx,ny,nz
    integer i,j,k
    real(KIND=CUSTOM_REAL) f(nx,ny,nz),f1(nx,ny,nz),h


    f1(:,:,:)=0.


    do k=2,nz-1
       do j=2,ny-1
          do i=2,nx-1

             f1(i,j,k) = f(i+1,j,k)-f(i-1,j,k)

          enddo
       enddo
    enddo


    f1(:,:,nz)=f1(:,:,nz-1)
    f1(:,ny,:)=f1(:,ny-1,:)
    f1(nx,:,:)=f1(nx-1,:,:)
    f1(:,:,1)=f1(:,:,2)
    f1(:,1,:)=f1(:,2,:)
    f1(1,:,:)=f1(2,:,:)

    f1(:,:,:) = f1(:,:,:) / (2.*h)


  end subroutine sub_D1


 subroutine sub_D2(f,f1,nx,ny,nz,h)

    implicit none
    integer nx,ny,nz
    integer i,j,k
    real(KIND=CUSTOM_REAL) h,f(nx,ny,nz),f1(nx,ny,nz)


    f1(:,:,:)=0.


    do k=2,nz-1
       do j=2,ny-1
          do i=2,nx-1

             f1(i,j,k) = f(i,j+1,k)-f(i,j-1,k)

          enddo
       enddo
    enddo


    f1(:,:,nz)=f1(:,:,nz-1)
    f1(:,ny,:)=f1(:,ny-1,:)
    f1(nx,:,:)=f1(nx-1,:,:)
    f1(:,:,1)=f1(:,:,2)
    f1(:,1,:)=f1(:,2,:)
    f1(1,:,:)=f1(2,:,:)

    f1(:,:,:) = f1(:,:,:) / (2.*h)


  end subroutine sub_D2

  subroutine sub_D3(f,f1,nx,ny,nz,h)

    implicit none
    integer nx,ny,nz
    integer i,j,k
    real(KIND=CUSTOM_REAL) f(nx,ny,nz),f1(nx,ny,nz),h


    f1(:,:,:)=0.


    do k=2,nz-1
       do j=2,ny-1
          do i=2,nx-1

             f1(i,j,k) = f(i,j,k+1)-f(i,j,k-1)

          enddo
       enddo
    enddo


    f1(:,:,nz)=f1(:,:,nz-1)
    f1(:,ny,:)=f1(:,ny-1,:)
    f1(nx,:,:)=f1(nx-1,:,:)
    f1(:,:,1)=f1(:,:,2)
    f1(:,1,:)=f1(:,2,:)
    f1(1,:,:)=f1(2,:,:)

    f1(:,:,:) = f1(:,:,:) / (2.*h)

  end subroutine sub_D3




  subroutine sub_D11(f,f1,nx,ny,nz,h)

    implicit none
    integer nx,ny,nz
    integer i,j,k
    real(KIND=CUSTOM_REAL) f(nx,ny,nz),f1(nx,ny,nz),h


    f1(:,:,:)=0.


    do k=2,nz-1
       do j=2,ny-1
          do i=2,nx-1

             f1(i,j,k) = f(i+1,j,k)+f(i-1,j,k)-2.*f(i,j,k)

          enddo
       enddo
    enddo


    f1(:,:,nz)=f1(:,:,nz-1)
    f1(:,ny,:)=f1(:,ny-1,:)
    f1(nx,:,:)=f1(nx-1,:,:)
    f1(:,:,1)=f1(:,:,2)
    f1(:,1,:)=f1(:,2,:)
    f1(1,:,:)=f1(2,:,:)

    f1(:,:,:) = f1(:,:,:) / (h**2)

  end subroutine sub_D11

  subroutine sub_D22(f,f1,nx,ny,nz,h)

    implicit none
    integer nx,ny,nz
    integer i,j,k
    real(KIND=CUSTOM_REAL) f(nx,ny,nz),f1(nx,ny,nz),h


    f1(:,:,:)=0.


    do k=2,nz-1
       do j=2,ny-1
          do i=2,nx-1

             f1(i,j,k) = f(i,j+1,k)+f(i,j-1,k)-2.*f(i,j,k)

          enddo
       enddo
    enddo


    f1(:,:,nz)=f1(:,:,nz-1)
    f1(:,ny,:)=f1(:,ny-1,:)
    f1(nx,:,:)=f1(nx-1,:,:)
    f1(:,:,1)=f1(:,:,2)
    f1(:,1,:)=f1(:,2,:)
    f1(1,:,:)=f1(2,:,:)

    f1(:,:,:) = f1(:,:,:) / (h**2)

  end subroutine sub_D22

  subroutine sub_D33(f,f1,nx,ny,nz,h)

    implicit none
    integer nx,ny,nz
    integer i,j,k
    real(KIND=CUSTOM_REAL) f(nx,ny,nz),f1(nx,ny,nz),h


    f1(:,:,:)=0.


    do k=2,nz-1
       do j=2,ny-1
          do i=2,nx-1

             f1(i,j,k) = f(i,j,k+1)+f(i,j,k-1)-2.*f(i,j,k)

          enddo
       enddo
       !write(*,*) k,f1(2,2,k)
    enddo


    f1(:,:,nz)=f1(:,:,nz-1)
    f1(:,ny,:)=f1(:,ny-1,:)
    f1(nx,:,:)=f1(nx-1,:,:)
    f1(:,:,1)=f1(:,:,2)
    f1(:,1,:)=f1(:,2,:)
    f1(1,:,:)=f1(2,:,:)

    f1(:,:,:) = f1(:,:,:) / (h**2)


  end subroutine sub_D33

  subroutine sub_D13(f,f1,nx,ny,nz,h1,h3)

    implicit none
    integer nx,ny,nz
    integer i,j,k
    real(KIND=CUSTOM_REAL) f(nx,ny,nz),f1(nx,ny,nz),h1,h3


    f1(:,:,:)=0.


    do k=2,nz-1
       do j=2,ny-1
          do i=2,nx-1

             f1(i,j,k) = f(i+1,j,k+1)+f(i-1,j,k-1)-f(i+1,j,k-1)-f(i-1,j,k+1)

          enddo
       enddo
    enddo


    f1(:,:,nz)=f1(:,:,nz-1)
    f1(:,ny,:)=f1(:,ny-1,:)
    f1(nx,:,:)=f1(nx-1,:,:)
    f1(:,:,1)=f1(:,:,2)
    f1(:,1,:)=f1(:,2,:)
    f1(1,:,:)=f1(2,:,:)

    f1(:,:,:) = f1(:,:,:) / (4.*h1*h3)

  end subroutine sub_D13


  subroutine sub_D12(f,f1,nx,ny,nz,h1,h2)

    implicit none
    integer nx,ny,nz
    integer i,j,k
    real(KIND=CUSTOM_REAL) f(nx,ny,nz),f1(nx,ny,nz),h1,h2


    f1(:,:,:)=0.


    do k=2,nz-1
       do j=2,ny-1
          do i=2,nx-1

             f1(i,j,k) = f(i+1,j+1,k)+f(i-1,j-1,k)-f(i+1,j-1,k)-f(i-1,j+1,k)

          enddo
       enddo
    enddo


    f1(:,:,nz)=f1(:,:,nz-1)
    f1(:,ny,:)=f1(:,ny-1,:)
    f1(nx,:,:)=f1(nx-1,:,:)
    f1(:,:,1)=f1(:,:,2)
    f1(:,1,:)=f1(:,2,:)
    f1(1,:,:)=f1(2,:,:)

    f1(:,:,:) = f1(:,:,:) / (4.*h1*h2)

  end subroutine sub_D12


    subroutine sub_D23(f,f1,nx,ny,nz,h2,h3)

    implicit none
    integer nx,ny,nz
    integer i,j,k
    real(KIND=CUSTOM_REAL) f(nx,ny,nz),f1(nx,ny,nz),h2,h3


    f1(:,:,:)=0.


    do k=2,nz-1
       do j=2,ny-1
          do i=2,nx-1

             f1(i,j,k) = f(i,j+1,k+1)+f(i,j-1,k-1)-f(i,j+1,k-1)-f(i,j-1,k+1)

          enddo
       enddo
    enddo


    f1(:,:,nz)=f1(:,:,nz-1)
    f1(:,ny,:)=f1(:,ny-1,:)
    f1(nx,:,:)=f1(nx-1,:,:)
    f1(:,:,1)=f1(:,:,2)
    f1(:,1,:)=f1(:,2,:)
    f1(1,:,:)=f1(2,:,:)

    f1(:,:,:) = f1(:,:,:) / (4.*h2*h3)

  end subroutine sub_D23






end module deriv_DF_mod
