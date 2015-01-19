subroutine dclisb(a, n, nud, n1, np, b, eps, dr, z, ier)

!  simultaneous linear equations with real symmetric positive definite *
!      band matrix by cholesky method.                                 *
!  parameters                                                          *
!    (1) a : 2-dim. array containing the matrix.                       *
!    (2) n : order of the matrix.                                      *
!    (3) nud : size of band's half width.                              *
!    (4) n1 : row size of the array a in the 'dimension' statement.    *
!    (5) b : 1-dim. array containing the right hand side vector.       *
!    (6) eps : parameter to check singurarity off the matrix           *
!              standard value = 1.0d-14                                *
!    (7) dr : 1-dim. working array.                                    *
!    (8) z : 1-dim. working array.                                     *
!    (9) ier : error code.                                             *
!  copy right   t. oguni   july 30 1989   version 1.0                  *
!***********************************************************************
  integer :: n, nud, n1, np, ier
  complex(kind(0d0)) :: a(n1,n), b(n), dr(n), z(n)
  real(kind(0d0)) :: eps

  complex(kind(0d0)) :: xx, s, sum, au, t
  real(kind(0d0)) :: eps1
  integer :: i ,m, j, k1, mj, i1, k
  !  check the input data
  ier = 0
  eps1 = 1.0d-14
  m = nud + 1
  if ((n <= 0) .or. (nud <= 0 ) .or. (n1 < m)) then
     ier = 2
     print *, '(subr. lisb) invalid argument. ', n, nud, n1
  endif
  if (eps <= 0.0) eps = eps1
  !  modified cholesky decomposition
  j = 1
  if (abs(a(m,1)) <= eps) then
     ier = 1
     print *, '(subr. lisb) singular at step # ', j
  endif
  dr(1) = dcmplx(1.0d0) / a(m,1)
  xx = a(m-1,2)
  a(m-1,2) = a(m-1,2) * dr(1)
  s = a(m,2) - xx * a(m-1,2)
  j = 2
  if (abs(s) <= eps) then
     ier = 1
     print *, '(subr. lisb) singular at step # ', j
  endif
  dr(2) = dcmplx(1.0d0) / s
  if (m < 3) then
     do j=3,n
        xx = a(1,j)
        a(1,j) = xx * dr(j-1)
        s = a(2,j) - xx * a(1,j)
        if (abs(s) <= eps) then
           ier = 1
           print *, ' (subr. lisb) singular at step # ', j
        endif
        dr(j) = dcmplx(1.0d0) / s
     enddo
  else
     do j=3,n
        k1 = 1
        if (j >= m) k1 = j - m + 1
        mj = m - j
        do i=k1+1,j-1
           sum = dcmplx(0.0d0)
           do k=k1,i-1
              sum = sum + a(m-i+k,i) * a(mj+k,j)
           enddo
           a(mj+i,j) = a(mj+i,j) - sum
        enddo
        sum = dcmplx(0.0d0)
        do i=k1,j-1
           xx = a(mj+i,j)
           au = xx * dr(i)
           sum = sum + xx *au
           a(mj+i,j) = au
        enddo
        t = a(m,j) - sum
        if (abs(t) <= eps) then
           ier = 1
           print *, ' (subr. lisb) singular at step # ', j
        endif
        dr(j) = dcmplx(1.0d0) / t
     enddo
  endif
  ! subtitution
  entry dcsbsub(a, n, nud, n1, np, b, eps, dr, z, ier)
  !  forward substitution
  m = nud + 1
  if (m < 3) then
     z(np) = b(np)
     do j=np+1,n
        z(j) = b(j) - a(1,j) * z(j-1)
     enddo
     b(n) = z(n) * dr(n)
  else
     z(np) = b(np)
     z(np+1) = b(np+1) - a(m-1,np+1) * z(np)
     do j=np+2,n
        if (j > np-1+m) then
           i1 = 1
        else
           i1 = np-1+m - j + 1
        endif
        sum = dcmplx(0.0d0)
        do  k=i1,m-1
           sum = sum + a(k,j) * z(j-m+k)
        enddo
        z(j) = b(j) - sum
     enddo
     do  j=n-1,n
        z(j) = z(j) * dr(j)
     enddo
     b(n) = z(n)
     b(n-1) = z(n-1) - a(m-1,n) * z(n)
  endif

end subroutine dclisb


