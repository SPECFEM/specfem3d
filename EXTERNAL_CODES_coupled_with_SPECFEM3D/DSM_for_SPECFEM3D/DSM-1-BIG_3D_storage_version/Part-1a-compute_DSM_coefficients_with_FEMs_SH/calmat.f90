
subroutine calmatc( nlayer,vnp,vra,con,rpow,w1dn,w2dn,ra,m,work )

  ! Computing \int r^rpow con W_p^(w1dn) W_q^(w2dn) dr.
  implicit none
  integer, parameter ::  maxrpow = 2
  integer :: nlayer,vnp,rpow,w1dn,w2dn
  real(kind(0d0)) :: vra(vnp),con(vnp),ra(*),m(*),work(*)
  integer :: i,j,k,kk,l,nn,snp
  real(kind(0d0)) :: a(2,2),b(2,2),c(5),rh
  ! parameter check
  if ( rpow > maxrpow ) pause 'Invalid arguments.(calmatc)'
  ! computing of the matrix elements
  snp = 1
  do i=1,nlayer
     rh = ra(i+1) - ra(i)
     if ( w1dn == 0 ) then
        a(2,1) = -1.d0 / rh
        a(1,1) = ra(i+1) / rh
        a(2,2) = 1.d0 / rh
        a(1,2) = -ra(i) / rh
     else
        if ( w1dn == 1 ) then
           a(2,1) = 0.d0
           a(1,1) = -1.d0 / rh
           a(2,2) = 0.d0
           a(1,2) = 1.d0 / rh
        else
           pause 'Invalid arguments.(calmatc)'
        endif
     endif
     if ( w2dn == 0 ) then
        b(2,1) = -1.d0 / rh
        b(1,1) = ra(i+1) / rh
        b(2,2) = 1.d0 / rh
        b(1,2) = -ra(i) / rh
     else
        if ( w2dn == 1 ) then
           b(2,1) = 0.d0
           b(1,1) = -1.d0 / rh
           b(2,2) = 0.d0
           b(1,2) = 1.d0 / rh
        else
           pause 'Invalid arguments.(calmatc)'
        endif
     endif
     do j=1,2
        do k=1,2
           do kk=1,5
              c(kk) = 0.d0
           enddo
           call pmulti( 2,a(1,j),2,b(1,k),3,c )
           do l=3,1,-1
              c(l+rpow) = c(l)
              if ( rpow > 0 ) c(l)=0.d0
           enddo
           nn = 4 * (i-1) + 2 * (j-1) + k
           call pinteg( snp,5,c,ra(i),ra(i+1), &
                &  vnp,vra,con,work(nn) )
        enddo
     enddo
  enddo

  if ( w1dn /= w2dn ) then
     do i=1,4*nlayer
        if ( (mod(i,4) == 0) .or. (mod(i,4) == 1) ) then
           m(i) = 2.d0 * work(i)
        else
           if ( mod(i,4) == 2 ) then
              m(i) = work(i) + work(i+1)
           else
              m(i) = work(i-1) + work(i)
           endif
        endif
     enddo
  else
     do i=1,4*nlayer
        m(i) = work(i)
     enddo
  endif


end subroutine calmatc

!-----------------------------------------------------------------------

subroutine caltl( nlayer,vnp,vra,rho,ra,tl)

  ! Computing of lumped mass matrix.
  implicit none
  integer :: nlayer,vnp
  real(kind(0d0)) :: vra(vnp),rho(vnp),ra(*),tl(*)
  integer :: i,nn,snp
  real(kind(0d0)) ::c(3),from,to

  snp = 1

  do i=1,nlayer
     c(1) = 0.d0
     c(2) = 0.d0
     c(3) = 1.d0

     from = ra(i)
     to = ( ra(i) + ra(i+1) ) / 2.d0
     nn = 4 * (i-1)
     call pinteg( snp,3,c,from,to,vnp,vra,rho,tl(nn+1) )

     tl(nn+2) = 0.d0
     tl(nn+3) = 0.d0

     from = to
     to = ra(i+1)
     call pinteg( snp,3,c,from,to,vnp,vra,rho,tl(nn+4) )

  enddo


end subroutine caltl

!-----------------------------------------------------------------------

subroutine calhl( nlayer,vnp,vra,mu,ra,hl)

  ! Computing of lumped rigidity matrix.
  implicit none
  integer :: nlayer,vnp
  real(kind(0d0)) :: vra(vnp),mu(vnp),ra(*),hl(*)
  integer ::i,nn,snp
  real(kind(0d0)) :: c(1),from,to

  snp = 1
  do i=1,nlayer
     c(1) = 1.d0

     from = ra(i)
     to = ( ra(i) + ra(i+1) ) / 2.d0
     nn = 4 * (i-1)
     call pinteg( snp,1,c,from,to,vnp,vra,mu,hl(nn+1) )

     hl(nn+2) = 0.d0
     hl(nn+3) = 0.d0

     from = to
     to = ra(i+1)
     call pinteg( snp,1,c,from,to,vnp,vra,mu,hl(nn+4) )
  enddo
end subroutine calhl

!-----------------------------------------------------------------------

subroutine calt( nlayer, tl, tc, t )
  implicit none
  integer ::  nlayer, i
  real(kind(0d0)) :: t(*), tl(*), tc(*)

  do i=1,4*nlayer
     t(i) = ( tl(i) + tc(i) ) / 2.d0
  enddo
end subroutine calt


!-----------------------------------------------------------------------

subroutine pmulti(n,a,m,b,l,c)

! Computing the (l-1) degrees polynomial c(n) which is the product of
! (n-1) degrees polynomial a(n) and (m-1) degrees polynomial b(n).
  implicit none
  integer :: n,m,l,i,j
  real(kind(0d0)) :: a(n),b(m),c(l)

  if (n+m-1 /= l) pause 'Invalid arguments.(pmulti)'
  do i=1,l
     c(i)=0.d0
  enddo

  do i=1,n
     do j=1,m
        c(i+j-1) = c(i+j-1) + a(i)*b(j)
     enddo
  enddo
end subroutine pmulti

!-----------------------------------------------------------------------

subroutine pinteg(snp,n,p,from,to,vnp,vra,con,pint)

! Evaluating the integrated value pint from 'from' to 'to' of p(n)*con
! which is the product of (n-1) degrees polynomial 'p(n)' and 'con'.
  implicit none
  integer, parameter :: maxn = 5
  integer :: snp,n,vnp
  real(kind(0d0)) :: from,to,p(n),vra(vnp),con(vnp),pint
  real(kind(0d0)) :: x1,x2,q(2),pq(maxn+1),psint


  if ( n > maxn ) pause 'Degrees of a polynomial is too large.(pinteg)'
  if ( snp >= vnp ) snp = 1

  pint = 0.d0
  x1 = from

100 continue
  if ( (vra(snp) <= x1) .and. (vra(snp+1) > x1) ) then
     x2 = dmin1( to, vra(snp+1) )
  else
     snp = snp + 1
     if ( snp == vnp ) snp = 1
     goto 100
  endif
  ! evaluating the integrated value
  q(2) = ( con(snp+1)-con(snp) ) / ( vra(snp+1)-vra(snp) )
  q(1) = - q(2) * vra(snp) + con(snp)
  call pmulti( n,p,2,q,n+1,pq )
  call polint( n+1,pq,x1,x2,psint )
  pint = pint + psint

  if ( x2 == to ) then
     continue
  else
     x1 = x2
     goto 100
  endif

end subroutine pinteg

!-----------------------------------------------------------------------


subroutine polint( n,p,x1,x2,pint )

  ! Evaluating the integrated value 'pint' from 'x1' to 'x2'
  ! of (n-1) degrees polynomial 'p(n)'.
  implicit none
  integer, parameter :: maxn=6
  integer :: n
  real(kind(0d0)) :: p(*),x1,x2,pint
  integer :: i,j
  real(kind(0d0)) :: a(maxn),b(maxn),dx,xx

  if ( n > maxn) pause 'Degrees of a polynomial is too large.(polint)'

  a(1) = 1.d0
  b(1) = 1.d0
  if ( n >= 2 ) then
     do i=2,n
        a(i) = a(i-1) * x1
        b(i) = b(i-1) * x2
     enddo
  endif
  dx = x2 - x1

  pint = 0.d0
  do i=1,n
     xx = 0.d0
     do j=1,i
        xx = xx + a(j) * b(i-j+1) / dble(i)
     enddo
     pint = pint + p(i) * dx * xx
  enddo

end subroutine polint


!-----------------------------------------------------------------------

subroutine cala0( nlayer,omega,omegai,t,h1,h2,h3,h4,coef,a0 )

  ! Computing the coefficient matrix 'a' in the solid part.
  implicit none
  integer :: nlayer
  real(kind(0d0)) :: omega,omegai
  real(kind(0d0)) :: t(*)
  real(kind(0d0)) :: h1(*),h2(*),h3(*),h4(*)
  complex(kind(0d0)) :: comega2,coef,a0(*)
  integer :: i
  real(kind(0d0)) :: h

  comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
  do i=1,4*nlayer
     h = h1(i) - h2(i) + h3(i) - 2.d0 * h4(i)
     a0(i) = comega2 * dcmplx( t(i) ) - coef * dcmplx( h )
  enddo

end subroutine cala0


!-----------------------------------------------------------------------

subroutine cala2( nlayer,h4,coef,a2 )

  ! Computing the coefficient matrix 'a' in the solid part.
  implicit none
  integer :: nlayer
  real(kind(0d0)) :: h4(*)
  complex(kind(0d0)) :: coef,a2(*)
  integer :: i

  do i=1,4*nlayer
     a2(i) = - coef * dcmplx( h4(i) )
  enddo


end subroutine cala2


!-----------------------------------------------------------------------

subroutine cala( nn,l,lda,a0,a2,a )

  ! Computing the coefficient matrix 'a' in the solid part.
  implicit none
  integer :: nn,l,lda
  complex(kind(0d0)) :: a0(lda,*),a2(lda,*),a(lda,*)
  integer :: i,j
  real(kind(0d0)) :: xl2

  xl2 = dble(l) * dble(l+1)
  do j=1,nn
     do i=1,2
        a(i,j) = a0(i,j) + dcmplx(xl2) * a2(i,j)
     enddo
  enddo
end subroutine cala


!-----------------------------------------------------------------------

subroutine calga(nlayer,omega,omegai,l,t,h1,h2,h3,h4,coef,a)

  ! Computing the coefficient matrix 'a' in the solid part.

  implicit none
  integer :: nlayer,l
  real(kind(0d0)) :: omega,omegai
  real(kind(0d0)) :: t(*)
  real(kind(0d0)) :: h1(*),h2(*),h3(*),h4(*)
  complex(kind(0d0)) :: comega2,coef,a(*)
  integer :: i
  real(kind(0d0)) :: h,xl2m2

  h = 0.d0
  do i=1,4*nlayer
     a(i) = 0.d0
  enddo
  comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
  xl2m2 = dble(l) * dble(l+1) -2.d0
  do i=1,4*nlayer
     h = h1(i) - h2(i) + h3(i) + xl2m2 * h4(i)
     a(i) = comega2 * dcmplx( t(i) ) - coef * dcmplx( h )
  enddo

end subroutine calga


!-----------------------------------------------------------------------


subroutine calga2(nlayer,omega,omegai,l,t,h1,h2,h3,coef,a)

  ! Computing the coefficient matrix 'a' in the solid part.

  implicit none
  integer :: nlayer,l
  real(kind(0d0)) :: omega,omegai
  real(kind(0d0)) :: t(*)
  real(kind(0d0)) :: h1(*),h2(*),h3(*)
  complex(kind(0d0)) :: comega2,coef,a(*)
  integer :: i
  real(kind(0d0)) :: h,xl2m1

  h = 0.d0
  do i=1,4*nlayer
     a(i) = 0.d0
  enddo
  comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
  xl2m1 = dble(l) * dble(l+1) -1.d0
  do i=1,4*nlayer
     h = h1(i) - h2(i) + xl2m1 * h3(i)
     a(i) = comega2 * dcmplx( t(i) ) - coef * dcmplx( h )
  enddo
end subroutine calga2

!-----------------------------------------------------------------------

subroutine overlap( nlayer,a,a2 )

  ! Overlapping the coefficient matrix elements in the solid part.
  implicit none
  integer :: nlayer
  complex(kind(0d0)) :: a(*),a2(2,*)
  integer :: i,j,k,mu,m,i1,i2,k1,k2,nsize


  nsize = nlayer+1
  mu = 1
  m = mu + 1

  do j=1,nsize
     i1 = max0(1,j-mu)
     i2 = j
     do i=i1,i2
        k = i - j + m
        if ( i == j ) then
           if ( i == 1 ) then
              a2(k,j) = a2(k,j) + a(1)
           else
              if ( i == nsize ) then
                 a2(k,j) = a2(k,j) + a( 4*nlayer )
              else
                 k1 = 4 * i - 4
                 k2 = k1 + 1
                 a2(k,j) = a2(k,j) + a(k1) + a(k2)
              endif
           endif
        endif
        if (i+1 == j) then
           k1 = 4 * i - 2
           a2(k,j) = a2(k,j) + a(k1)
        endif
     enddo
  enddo
end subroutine overlap

!-----------------------------------------------------------------------

subroutine calg2( l,m,spo,r0,mt,mu0,coef,ga,a,ga2,dr,g2,ig2 )

  implicit none
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  integer :: l,m
  real(kind(0d0)) :: spo,r0,mt(3,3),mu0
  complex(kind(0d0)) :: coef,ga(8),a(4),ga2(2,3),g2(*)
  integer :: i,itmp,ig2
  real(kind(0d0)) :: b1,sgn,ier
  complex(kind(0d0)) :: dd,cg2(3),dr(3),z(3)
  real(kind(0d0)) :: eps
  data eps / -1.d0 /


  ! computing of particular solution with free surface boundary conditions
  cg2 = 0
  dd = dcmplx( 0.d0 )
  if ( m >= 0 ) then
     sgn = 1.d0
  else
     sgn = - 1.d0
  endif
  if ( iabs(m) == 1 ) then
     b1 = dsqrt( dble( 2*l+1 ) / ( 16.d0 * pi ) )

     dd = dcmplx( b1 ) * ( dcmplx( sgn * mt(1,3), mt(1,2) ) )&
          &        / ( dcmplx( r0 * r0 * mu0 ) * coef )
     itmp = 4
     cg2 = cmplx(0.d0)
     do i=2,3
        cg2(i) = - dd * ( ga(itmp+1) + ga(itmp+2) )
        itmp = itmp + 2
     enddo
  else
     if ( iabs(m) == 2 ) then
        b1 = dsqrt( dble(2*l+1)*dble(l-1)*dble(l+2)/(64.d0*pi) )
        cg2(2) = dcmplx( b1 / r0 ) &
             &  * dcmplx( 2.d0*mt(2,3), sgn*( mt(2,2)-mt(3,3) ) )
     endif
  endif
  if (  ((m == -2) .or. (m == -l)) .and. (ig2 == 0) ) then
     call dclisb0( ga2,3,1,2,cg2,eps,dr,z,ier)
     ig2= ig2+1
  else
     call dcsbsub0( ga2,3,1,2,cg2,eps,dr,z,ier)
  endif
  cg2(3) = cg2(3) + dd
  ! computation of the excitation vector
  itmp = dint(spo)
  g2(itmp+1) = a(1) * cg2(1) + a(2) * cg2(3)
  g2(itmp+2) = a(3) * cg2(1) + a(4) * cg2(3)

end subroutine calg2
