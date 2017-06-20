
subroutine glu( a, n, n1, b, eps, wk, ip, ier )
  !
  !        glu, glusub
  !               copyright : h.hasegawa, aug. 26 1989 v.1
  !
  !               solves simultaneous linear equations
  !               by Gaussian elimination method.
  !
  !        input - -
  !             a(n1,n)  r *8  : 2-dim. array containing the coefficients.
  !             n        i *4  : order of matrix.
  !             n1       i *4  : size of array a.
  !             b(n)     r *8  : 1-dim. array containing the right-hand
  !                              side vector.
  !             eps      r *8  : parameter to check singularity of the
  !                              matrix. ( standard value 3.52d-15 )
  !        output - -
  !             a(n1,n)        : result of Gaussian elimination.
  !             b(n)           : solution.
  !             ip(n)    i *4  : pivot number.
  !             ier      i *4  : = 0,  for normal execution.
  !                              = 1,  for singular matrix.
  !                              = 2,  for singular original matrix.
  !                              = 3,  for invalid arguement.
  !        working  -
  !             wk(n)    r *8  : 1-dim. array.

  implicit none
  integer:: n,n1
  integer:: ip(n),ier
  real(kind(0d0)):: eps
  complex(kind(0d0)):: a(n1,n),b(n),wk(n)
  integer:: i,j,k,ipk
  complex(kind(0d0)):: amax,aik,w,t
  !             left-hand side
  if ( eps < 0.0d0 )  eps = 3.52d-15
  if ( ( n1 < n ) .or. ( n <= 0 ) ) then
     ier = 3
     print *, '  (subr. glu)  invalid argument.  n1, n =', n1, n
     return
  endif
  !             check original matrix.
  do i = 1, n
     wk(i) = zabs(a(i,1))
  enddo
  do j = 2, n
     do i = 1, n
        wk(i) = dmax1( zabs( wk(i) ), zabs(a(i,j)) )
     enddo
  enddo
  do i = 1, n
     if ( zabs( wk(i) ) < eps ) then
        ier = 2
        print *, '  (subr. glu)  original matrix is singular.'
        return
     endif
  enddo

  ier = 0
  do k = 1, n
     !             find maximum element in the k-th column.
     amax = zabs(a(k,k))
     ipk = k
     do i = k+1, n
        aik = zabs(a(i,k))
        if ( zabs( aik ) > zabs( amax ) ) then
           ipk = i
           amax = aik
        endif
     enddo
     ip(k) = ipk
     !
     if ( zabs( amax ) > eps ) then
        if ( ipk /= k ) then
           w = a(ipk,k)
           a(ipk,k) = a(k,k)
           a(k,k) = w
        endif
        !             compute alfa
        do i = k+1, n
           a(i,k) = -a(i,k)/a(k,k)
           wk(i) = a(i,k)
        enddo

        do j = k+1, n
           if ( ipk /= k ) then
              w = a(ipk,j)
              a(ipk,j) = a(k,j)
              a(k,j) = w
           endif
           !             Gaussian elimination
           t = a(k,j)
           do i = k+1, n
              a(i,j) = a(i,j) + wk(i)*t
           enddo
        enddo

        !             matrix is singular.
     else
        ier = 1
        ip(k) = k
        do i = k+1, n
           a(i,k) = 0.0d0
        enddo
        print *, '  (subr. glu)  matrix is singular at k =', k
        return
     endif
  enddo
  !             right-hand side
  ! entry glusub( a, b )
  ! entry glusub( a, n, n1, b, eps, wk, ip, ier )
  !           forward elimination process
  do k = 1, n
     if ( ip(k) /= k ) then
        w = b(ip(k))
        b(ip(k)) = b(k)
        b(k) = w
     endif

     t = b(k)
     do i = k+1, n
        b(i) = b(i) + a(i,k)*t
     enddo
  enddo

  !             backward substitution process
  b(n) = b(n)/a(n,n)
  do k = n-1, 1, -1
     t = b(k+1)
     do i = 1, k
        b(i) = b(i) - a(i,k+1)*t
     enddo
     b(k) = b(k)/a(k,k)
  enddo
  return
end subroutine glu

!---------------------------------------------------


subroutine interpolate( ncomp,nderiv,rsta,rrsta,g,u )

  implicit none
  integer:: ncomp,nderiv
  real(kind(0d0)):: rsta,rrsta(3)
  complex(kind(0d0)):: g(3*ncomp),u(ncomp)
  real(kind(0d0)):: dh(3)

  integer:: ip(3),ier,i,itmp,icomp
  complex(kind(0d0)):: a(3,3),b(3),wk(3)
  real(kind(0d0)):: eps
  eps = -1.d0

  do icomp=1,ncomp
     u(icomp) = dcmplx(0.d0)
  enddo
  do i=1,3
     dh(i) = rrsta(i) - rsta
  enddo

  if ( (dh(2) == 0.d0) .and. (nderiv == 0)) then
     itmp = ncomp + 1
     do icomp=1,ncomp
        u(icomp) = g(itmp)
        itmp = itmp + 1
     enddo
     return
  endif
  do i=1,3
     a(1,i) = dcmplx( 1.d0 )
     a(2,i) = dcmplx( dh(i) )
     a(3,i) = dcmplx( dh(i) * dh(i) / 2.d0 )
  enddo
  call fillinpb(nderiv,b)
  call glu(a,3,3,b,eps,wk,ip,ier)

  do icomp=1,ncomp
     do i=1,3
        u(icomp) = u(icomp) + b(i) * g( ncomp * (i-1) + icomp )
     enddo
  enddo

end subroutine interpolate

