subroutine rk3(neq,func,x0,xe,n,y0,yn,n1,work)
  !**********************************************************************
  !     subroutine rk numerically integrates a system of neq           *
  !     first order ordinary differential equations of the form        *
  !             dy(i)/dx = f(x, y(1),..., y(neq)),                     *
  !     by the classical runge-kutta formula.                          *
  !                                                                    *
  !     parameters                                                     *
  !  === input ===                                                     *
  !     (1) neq: number of equations to be integrated                  *
  !     (2) func: subroutine func(x,y,f) to evaluate derivatives       *
  !                f(i)=dy(i)/dx                                       *
  !     (3) x0: initial value of independent variable                  *
  !     (4) xe: output point at which the solution is desired          *
  !     (5) n: number of divisions                                     *
  !        the interval (x0, xe) is divided into n subintervals        *
  !        with the length (xe-x0)/n and in each subinterval           *
  !        the classical runge-kutta formula is used.                  *
  !     (6) y0(i) (i=1,..,neq): initial value at x0                    *
  !  === output ===                                                    *
  !     (7) yn(i) (i=1,..,neq): approximate solution at xe             *
  !  === other ===                                                     *
  !     (8) work(): two-dimentional array (size=(neq,2)) to be         *
  !                 used inside rk                                     *
  !     copyright: m. sugihara, november 15, 1989, v. 1                *
  !*********************************************************************
  implicit none
  external func

  integer:: neq,n,i,j,n1
  real(kind(0d0)):: x0,xe
  complex(kind(0d0)):: y0(neq),yn(neq),work(n1,2)
  real(kind(0d0)):: h
  h = (xe - x0) / dble(n)
  do i = 1,n
     call rkstep(neq,func,x0,h,y0,yn,work(1,1),work(1,2))
     x0 = x0 + h
     do j = 1,neq
        y0(j) = yn(j)
     enddo
  enddo
  x0 = xe
  return
end subroutine rk3

subroutine rkstep(neq,func,x,h,y0,yn,ak,w)
  implicit none
  real(kind(0d0)),parameter:: a2 =0.5d0,a3=a2,b2=0.5d0,b3=b2,c1=1.d0/6.d0,c2=1.d0/3.d0,c3=c2,c4=c1
  ! parameter(a2 = 0.5d0, a3 = a2)
  ! parameter(b2 = 0.5d0, b3 = b2)
  ! parameter(c1 = 1.d0/6.d0, c2 = 1.d0/3.d0, c3 = c2, c4 = c1)
  integer:: neq,i
  real(kind(0d0)):: x,h
  complex(kind(0d0)):: y0(neq),yn(neq),ak(neq),w(neq)
  external func
  call func(x,y0,ak)
  do i = 1,neq
     yn(i) = y0(i) + dcmplx( h * c1 ) * ak(i)
  enddo
  do i = 1,neq
     w(i) = y0(i) + dcmplx( h * b2 ) * ak(i)
  enddo
  call func(x + a2 * h,w,ak)
  do i = 1,neq
     yn(i) = yn(i) + dcmplx( h * c2 ) * ak(i)
  enddo
  do i = 1,neq
     w(i) = y0(i) + dcmplx( h * b3 ) * ak(i)
  enddo
  call func(x + a3 * h,w,ak)
  do i = 1,neq
     yn(i) = yn(i) + dcmplx( h * c3 ) * ak(i)
  enddo
  do i = 1,neq
     w(i) = y0(i) + dcmplx( h ) * ak(i)
  enddo
  call func(x + h,w,ak)
  do i = 1,neq
     yn(i) = yn(i) + dcmplx( h * c4 ) * ak(i)
  enddo
  return

end subroutine rkstep
