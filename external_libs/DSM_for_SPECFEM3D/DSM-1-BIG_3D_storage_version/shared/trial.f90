
!! DK DK when l is greater than 5 we can get rid of all the "if" tests in this subroutine to make it much faster

subroutine caldvec_for_l_more_than_5_no_store_l(l,theta,phi,plm,dvec,dvecdt,dvecdp,theta_n)

!! DK DK for Vadim: in order to avoid array memory copies, I have now put the loop on 'itheta' inside this subroutine.
!! DK DK for Vadim: I thus slightly changed its name to avoid confusion, because the list of arguments has changed.

!! DK DK for Vadim: for array plm() I also move index 'itheta' as first index; thus please do it also in all the rest of the code

  implicit none

  real(kind(0d0)),parameter :: pi=3.1415926535897932d0
  real(kind(0d0)),parameter :: convert_to_radians = pi/180.d0

  integer :: l,m,i,theta_n,itheta

  real(kind(0d0)):: theta(theta_n),phi(theta_n),plm(theta_n,1:3,0:3)
  real(kind(0d0)):: x,fact,thetaval,phival,inv_sintheta
  complex(kind(0d0)) :: dvec(theta_n,1:3,-2:2),expimp
  complex(kind(0d0)) :: dvecdt(theta_n,1:3,-2:2),dvecdp(theta_n,1:3,-2:2)
  real(kind(0d0)):: plmdt,xl2

!! DK DK modified this to precompute the coefficients
  real(kind(0d0)):: coef(0:2)

!! DK DK added this to precompute the coefficients
! do m=0,min(l,2)
  do m=0,2
     fact = 1.d0
     if ( m /= 0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef(m) = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
  enddo

!! DK DK added this loop (it was previously in the calling program)
!! DK DK adds an Intel compiler directive to force vectorization here, otherwise the Intel compiler
!! DK DK displays this for some reason: "remark: loop was not vectorized: vectorization possible but seems inefficient".
!! DK DK Also added a compiler directive for the xlf compiler on IBM Blue Gene
!! from http://pic.dhe.ibm.com/infocenter/compbg/v121v141/index.jsp?topic=%2Fcom.ibm.xlf141.bg.doc%2Flanguage_ref%2Fassert.html
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
  do itheta = 1,theta_n

!! DK DK added this
  thetaval = theta(itheta) * convert_to_radians
  phival = phi(itheta) * convert_to_radians

!! DK DK added this
!! DK DK it is always better to precompute an inverse and store it if we use it several times in the loop,
!! DK DK because on processors divisions are much more expensive than multiplications and thus it is much better
!! DK DK to compute the inverse once and for all and then later multiply by the inverse that has been stored
  inv_sintheta = 1.d0 / dsin(thetaval)

  x = dcos( thetaval )
  xl2 = dble(l) * dble(l+1)

!! DK DK l starts at 0, thus the upper bound of the loop can be 0, 1, 2, or 3
! do m=0,min(l,3)
!   call calplm( l,m,x,plm(itheta,1,m) )
!! DK DK unrolled the loop, otherwise it prevents vectorization of the big outer loop on itheta.
!! DK DK also inlined the call to the subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 0 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     plm(itheta,3,0) = plm(itheta,2,0)
     plm(itheta,2,0) = plm(itheta,1,0)
     plm(itheta,1,0) = (x*dble(2*l-1) * plm(itheta,2,0)-dble(l-1) * plm(itheta,3,0) )/dble(l)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 1 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     plm(itheta,3,1) = plm(itheta,2,1)
     plm(itheta,2,1) = plm(itheta,1,1)
     plm(itheta,1,1) = (x*dble(2*l-1) * plm(itheta,2,1)-dble(l) * plm(itheta,3,1) )/dble(l-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 2 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     plm(itheta,3,2) = plm(itheta,2,2)
     plm(itheta,2,2) = plm(itheta,1,2)
     plm(itheta,1,2) = (x*dble(2*l-1) * plm(itheta,2,2)-dble(l+1) * plm(itheta,3,2) )/dble(l-2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 3 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     plm(itheta,3,3) = plm(itheta,2,3)
     plm(itheta,2,3) = plm(itheta,1,3)
     plm(itheta,1,3) = (x*dble(2*l-1) * plm(itheta,2,3)-dble(l+2) * plm(itheta,3,3) )/dble(l-3)

! enddo


!! DK DK l starts at 0, thus the upper bound of the loop can be 0, 1, or 2
! do m=0,min(l,2)
!! DK DK unrolled the loop, otherwise it prevents vectorization of the big outer loop on itheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 0 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    expimp = zexp( dcmplx( 0.d0, 0.d0*phival ) )
     expimp = zexp( dcmplx( 0.d0, 0.d0 ) )
!    plmdt = 0.d0 * x * inv_sintheta * plm(itheta,1,0) + plm(itheta,1,1)
     plmdt = plm(itheta,1,1)
     dvec(itheta,1,0)  = coef(0) * plm(itheta,1,0) * expimp
!    dvec(itheta,1,-0) = dconjg( dvec(itheta,1,0) )
     dvec(itheta,2,0) = coef(0) * plmdt * expimp
!    dvec(itheta,2,-0) = dconjg( dvec(itheta,2,0) )
     dvec(itheta,3,0)  = dcmplx( 0.d0, 0.d0) !!!! *inv_sintheta*coef(0)*plm(itheta,1,0)*expimp
!    dvec(itheta,3,-0) = dconjg( dvec(itheta,3,0))

     ! calculate derivatives
     dvecdt(itheta,1,0) = plmdt * coef(0) * expimp
!    dvecdt(itheta,1,-0) = dconjg( dvecdt(itheta,1,0) )
!    dvecdt(itheta,2,0) =(-x*inv_sintheta*plmdt + 0.d0*0.d0/(1-x*x)*plm(itheta,1,0) - xl2*plm(itheta,1,0)) *coef(0)*expimp
     dvecdt(itheta,2,0) =(-x*inv_sintheta*plmdt - xl2*plm(itheta,1,0)) *coef(0)*expimp
!    dvecdt(itheta,2,-0) = dconjg( dvecdt(itheta,2,0) )
     dvecdt(itheta,3,0) = dcmplx( 0.d0, 0.d0 ) !!! * ( - x / ( 1- x * x ) *plm(itheta,1,0) + inv_sintheta*plmdt)*coef(0)*expimp
!    dvecdt(itheta,3,-0) = dconjg( dvecdt(itheta,3,0) )
     dvecdp(itheta,1,0) = dcmplx( 0.d0, 0.d0 ) !!! * plm(itheta,1,0) * coef(0) * expimp
!    dvecdp(itheta,1,-0) = dconjg( dvecdp(itheta,1,0) )
     dvecdp(itheta,2,0) = dcmplx( 0.d0,0.d0 ) !!! * plmdt * coef(0) * expimp
!    dvecdp(itheta,2,-0) = dconjg( dvecdp(itheta,2,0) )
!    dvecdp(itheta,3,0) = - 0.d0 * 0.d0 *inv_sintheta*plm(itheta,1,0) * coef(0) * expimp
     dvecdp(itheta,3,0) = dcmplx( 0.d0,0.d0 )
!    dvecdp(itheta,3,-0) = dconjg( dvecdp(itheta,3,0) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 1 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     expimp = zexp( dcmplx( 0.d0, phival ) )
     plmdt = x * inv_sintheta * plm(itheta,1,1) + plm(itheta,1,2)
     dvec(itheta,1,1)  = coef(1) * plm(itheta,1,1) * expimp
     dvec(itheta,1,-1) = - dconjg( dvec(itheta,1,1) )
     dvec(itheta,2,1) = coef(1) * plmdt * expimp
     dvec(itheta,2,-1) = - dconjg( dvec(itheta,2,1) )
     dvec(itheta,3,1)  = dcmplx( 0.d0, 1.d0)*inv_sintheta*coef(1)*plm(itheta,1,1)*expimp
     dvec(itheta,3,-1) = - dconjg( dvec(itheta,3,1))

     ! calculate derivatives
     dvecdt(itheta,1,1) = plmdt * coef(1) * expimp
     dvecdt(itheta,1,-1) = - dconjg( dvecdt(itheta,1,1) )
     dvecdt(itheta,2,1) =(-x*inv_sintheta*plmdt+1.d0/(1-x*x)*plm(itheta,1,1)-xl2*plm(itheta,1,1))*coef(1)*expimp
     dvecdt(itheta,2,-1) = - dconjg( dvecdt(itheta,2,1) )
     dvecdt(itheta,3,1) = dcmplx( 0.d0, 1.d0 ) * ( - x / ( 1- x * x ) *plm(itheta,1,1) + inv_sintheta*plmdt)*coef(1)*expimp
     dvecdt(itheta,3,-1) = - dconjg( dvecdt(itheta,3,1) )
     dvecdp(itheta,1,1) = dcmplx( 0.d0, 1.d0 ) * plm(itheta,1,1) * coef(1) * expimp
     dvecdp(itheta,1,-1) = - dconjg( dvecdp(itheta,1,1) )
     dvecdp(itheta,2,1) = dcmplx( 0.d0,1.d0 ) * plmdt * coef(1) * expimp
     dvecdp(itheta,2,-1) = - dconjg( dvecdp(itheta,2,1) )
     dvecdp(itheta,3,1) = - inv_sintheta*plm(itheta,1,1) * coef(1) * expimp
     dvecdp(itheta,3,-1) = - dconjg( dvecdp(itheta,3,1) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 2 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     expimp = zexp( dcmplx( 0.d0, 2.d0*phival ) )
     plmdt = 2.d0 * x * inv_sintheta * plm(itheta,1,2) + plm(itheta,1,3)
     dvec(itheta,1,2)  = coef(2) * plm(itheta,1,2) * expimp
     dvec(itheta,1,-2) = dconjg( dvec(itheta,1,2) )
     dvec(itheta,2,2) = coef(2) * plmdt * expimp
     dvec(itheta,2,-2) = dconjg( dvec(itheta,2,2) )
     dvec(itheta,3,2)  = dcmplx( 0.d0, 2.d0)*inv_sintheta*coef(2)*plm(itheta,1,2)*expimp
     dvec(itheta,3,-2) = dconjg( dvec(itheta,3,2))

     ! calculate derivatives
     dvecdt(itheta,1,2) = plmdt * coef(2) * expimp
     dvecdt(itheta,1,-2) = dconjg( dvecdt(itheta,1,2) )
     dvecdt(itheta,2,2) =(-x*inv_sintheta*plmdt+4.d0/(1-x*x)*plm(itheta,1,2)-xl2*plm(itheta,1,2))*coef(2)*expimp
     dvecdt(itheta,2,-2) = dconjg( dvecdt(itheta,2,2) )
     dvecdt(itheta,3,2) = dcmplx( 0.d0, 2.d0 ) * ( - x / ( 1- x * x ) *plm(itheta,1,2) + inv_sintheta*plmdt)*coef(2)*expimp
     dvecdt(itheta,3,-2) = dconjg( dvecdt(itheta,3,2) )
     dvecdp(itheta,1,2) = dcmplx( 0.d0, 2.d0 ) * plm(itheta,1,2) * coef(2) * expimp
     dvecdp(itheta,1,-2) = dconjg( dvecdp(itheta,1,2) )
     dvecdp(itheta,2,2) = dcmplx( 0.d0,2.d0 ) * plmdt * coef(2) * expimp
     dvecdp(itheta,2,-2) = dconjg( dvecdp(itheta,2,2) )
     dvecdp(itheta,3,2) = - 4.d0 *inv_sintheta*plm(itheta,1,2) * coef(2) * expimp
     dvecdp(itheta,3,-2) = dconjg( dvecdp(itheta,3,2) )

! enddo

!! DK DK added this loop (it was previously in the calling program)
  enddo

end subroutine caldvec_for_l_more_than_5_no_store_l

!---------------------------------------------------------------

subroutine caldvec_for_l_less_than_5_no_store_l(l,theta,phi,plm,dvec,dvecdt,dvecdp,theta_n)

!! DK DK for Vadim: in order to avoid array memory copies, I have now put the loop on 'itheta' inside this subroutine.
!! DK DK for Vadim: I thus slightly changed its name to avoid confusion, because the list of arguments has changed.

!! DK DK for Vadim: for array plm() I also move index 'itheta' as first index; thus please do it also in all the rest of the code

  implicit none

  real(kind(0d0)),parameter :: pi=3.1415926535897932d0
  real(kind(0d0)),parameter :: convert_to_radians = pi/180.d0

  integer :: l,m,i,theta_n,itheta

  real(kind(0d0)):: theta(theta_n),phi(theta_n),plm(theta_n,1:3,0:3)
  real(kind(0d0)):: x,fact,thetaval,phival,inv_sintheta
  complex(kind(0d0)) :: dvec(theta_n,1:3,-2:2),expimp
  complex(kind(0d0)) :: dvecdt(theta_n,1:3,-2:2),dvecdp(theta_n,1:3,-2:2)
  real(kind(0d0)):: plmdt,xl2

!! DK DK modified this to precompute the coefficients
  real(kind(0d0)):: coef(0:2)

!! DK DK added this to inline the calplm() routine
! real(kind(0d0)) :: pmm

!! DK DK added this to precompute the coefficients
  do m=0,min(l,2)
     fact = 1.d0
     if ( m /= 0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef(m) = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
  enddo

!! DK DK added this loop (it was previously in the calling program)
!! DK DK adds an Intel compiler directive to force vectorization here, otherwise the Intel compiler
!! DK DK displays this for some reason: "remark: loop was not vectorized: vectorization possible but seems inefficient".
!! DK DK Also added a compiler directive for the xlf compiler on IBM Blue Gene
!! from http://pic.dhe.ibm.com/infocenter/compbg/v121v141/index.jsp?topic=%2Fcom.ibm.xlf141.bg.doc%2Flanguage_ref%2Fassert.html
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
  do itheta = 1,theta_n

!! DK DK added this
  thetaval = theta(itheta) * convert_to_radians
  phival = phi(itheta) * convert_to_radians

!! DK DK added this
!! DK DK it is always better to precompute an inverse and store it if we use it several times in the loop,
!! DK DK because on processors divisions are much more expensive than multiplications and thus it is much better
!! DK DK to compute the inverse once and for all and then later multiply by the inverse that has been stored
  inv_sintheta = 1.d0 / dsin(thetaval)

  x = dcos( thetaval )
  xl2 = dble(l) * dble(l+1)

!! DK DK l starts at 0, thus the upper bound of the loop can be 0, 1, 2, or 3
! do m=0,min(l,3)
!   call calplm( l,m,x,plm(itheta,1,m) )
!! DK DK unrolled the loop, otherwise it prevents vectorization of the big outer loop on itheta.
!! DK DK also inlined the call to the subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 0 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ( l == 0 ) then
!    pmm = 1.d0
!! DK DK useless because m == 0 here
!    if ( m>0 ) then
!       somx2 = dsqrt( (1.d0-x)*(1.d0+x) )
!       fact = 1.d0
!       do i=1,m
!          pmm = -pmm * fact * somx2
!          fact = fact + 2.d0
!       enddo
!    endif
     plm(itheta,3,0) = 0.d0
     plm(itheta,2,0) = 0.d0
     plm(itheta,1,0) = 1.d0 !! DK DK   pmm
  else
     plm(itheta,3,0) = plm(itheta,2,0)
     plm(itheta,2,0) = plm(itheta,1,0)
     if ( l == 1 ) then
        plm(itheta,1,0) = x * plm(itheta,2,0)
     else
        plm(itheta,1,0) = (x*dble(2*l-1) * plm(itheta,2,0)-dble(l-1) * plm(itheta,3,0) )/dble(l)
     endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 1 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (l >= 1) then

  if ( l == 1 ) then
!    if ( 1>0 ) then
!       somx2 = dsqrt( (1.d0-x)*(1.d0+x) )
!       fact = 1.d0
!! DK DK avoid a small loop here because it would prevent vectorization of the big loop on itheta
!       do i=1,1
!          pmm = - dsqrt( (1.d0-x)*(1.d0+x) )   !! DK DK computed this analytically instead of with a loop
!          fact = fact + 2.d0
!       enddo
!    else
!      pmm = 1.d0
!    endif
     plm(itheta,3,1) = 0.d0
     plm(itheta,2,1) = 0.d0
     plm(itheta,1,1) = - dsqrt( (1.d0-x)*(1.d0+x) )   !! DK DK computed this analytically instead of with a loop
  else
     plm(itheta,3,1) = plm(itheta,2,1)
     plm(itheta,2,1) = plm(itheta,1,1)
     if ( l == 2 ) then
        plm(itheta,1,1) = x * 3.d0 * plm(itheta,2,1)
     else
        plm(itheta,1,1) = (x*dble(2*l-1) * plm(itheta,2,1)-dble(l) * plm(itheta,3,1) )/dble(l-1)
     endif
  endif

  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 2 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (l >= 2) then

  if ( l == 2 ) then
!    if ( 2>0 ) then
!       somx2 = dsqrt( (1.d0-x)*(1.d0+x) )
!       fact = 1.d0
!! DK DK avoid a small loop here because it would prevent vectorization of the big loop on itheta
!       do i=1,2
!          fact = fact + 2.d0
!          pmm = 3.d0 * dsqrt( (1.d0-x)*(1.d0+x) )**2   !! DK DK computed this analytically instead of with a loop
!          fact = fact + 2.d0
!       enddo
!    else
!      pmm = 1.d0
!    endif
     plm(itheta,3,2) = 0.d0
     plm(itheta,2,2) = 0.d0
     plm(itheta,1,2) = 3.d0 * dsqrt( (1.d0-x)*(1.d0+x) )**2   !! DK DK computed this analytically instead of with a loop
  else
     plm(itheta,3,2) = plm(itheta,2,2)
     plm(itheta,2,2) = plm(itheta,1,2)
     if ( l == 3 ) then
        plm(itheta,1,2) = x * 5.d0 * plm(itheta,2,2)
     else
        plm(itheta,1,2) = (x*dble(2*l-1) * plm(itheta,2,2)-dble(l+1) * plm(itheta,3,2) )/dble(l-2)
     endif
  endif

  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 3 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (l >= 3) then

  if ( l == 3 ) then
!    if ( 3>0 ) then
!       somx2 = dsqrt( (1.d0-x)*(1.d0+x) )
!       fact = 1.d0
!! DK DK avoid a small loop here because it would prevent vectorization of the big loop on itheta
!       do i=1,3
!          pmm = - 15.d0 * dsqrt( (1.d0-x)*(1.d0+x) )**3   !! DK DK computed this analytically instead of with a loop
!          fact = fact + 2.d0
!          pmm = -pmm * fact * somx2
!          fact = fact + 2.d0
!          pmm = -pmm * fact * somx2
!          fact = fact + 2.d0
!       enddo
!    else
!      pmm = 1.d0
!    endif
     plm(itheta,3,3) = 0.d0
     plm(itheta,2,3) = 0.d0
     plm(itheta,1,3) = - 15.d0 * dsqrt( (1.d0-x)*(1.d0+x) )**3   !! DK DK computed this analytically instead of with a loop
  else
     plm(itheta,3,3) = plm(itheta,2,3)
     plm(itheta,2,3) = plm(itheta,1,3)
     if ( l == 4 ) then
        plm(itheta,1,3) = x * 7.d0 * plm(itheta,2,3)
     else
        plm(itheta,1,3) = (x*dble(2*l-1) * plm(itheta,2,3)-dble(l+2) * plm(itheta,3,3) )/dble(l-3)
     endif
  endif

  endif

! enddo



!! DK DK l starts at 0, thus the upper bound of the loop can be 0, 1, or 2
! do m=0,min(l,2)
!! DK DK unrolled the loop, otherwise it prevents vectorization of the big outer loop on itheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 0 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    expimp = zexp( dcmplx( 0.d0, 0.d0*phival ) )
     expimp = zexp( dcmplx( 0.d0, 0.d0 ) )
!    plmdt = 0.d0 * x * inv_sintheta * plm(itheta,1,0) + plm(itheta,1,1)
     plmdt = plm(itheta,1,1)
     dvec(itheta,1,0)  = coef(0) * plm(itheta,1,0) * expimp
!    dvec(itheta,1,-0) = dconjg( dvec(itheta,1,0) )
     dvec(itheta,2,0) = coef(0) * plmdt * expimp
!    dvec(itheta,2,-0) = dconjg( dvec(itheta,2,0) )
     dvec(itheta,3,0)  = dcmplx( 0.d0, 0.d0) !!!! *inv_sintheta*coef(0)*plm(itheta,1,0)*expimp
!    dvec(itheta,3,-0) = dconjg( dvec(itheta,3,0))

     ! calculate derivatives
     dvecdt(itheta,1,0) = plmdt * coef(0) * expimp
!    dvecdt(itheta,1,-0) = dconjg( dvecdt(itheta,1,0) )
!    dvecdt(itheta,2,0) =(-x*inv_sintheta*plmdt + 0.d0*0.d0/(1-x*x)*plm(itheta,1,0) - xl2*plm(itheta,1,0)) *coef(0)*expimp
     dvecdt(itheta,2,0) =(-x*inv_sintheta*plmdt - xl2*plm(itheta,1,0)) *coef(0)*expimp
!    dvecdt(itheta,2,-0) = dconjg( dvecdt(itheta,2,0) )
     dvecdt(itheta,3,0) = dcmplx( 0.d0, 0.d0 ) !!! * ( - x / ( 1- x * x ) *plm(itheta,1,0) + inv_sintheta*plmdt)*coef(0)*expimp
!    dvecdt(itheta,3,-0) = dconjg( dvecdt(itheta,3,0) )
     dvecdp(itheta,1,0) = dcmplx( 0.d0, 0.d0 ) !!! * plm(itheta,1,0) * coef(0) * expimp
!    dvecdp(itheta,1,-0) = dconjg( dvecdp(itheta,1,0) )
     dvecdp(itheta,2,0) = dcmplx( 0.d0,0.d0 ) !!! * plmdt * coef(0) * expimp
!    dvecdp(itheta,2,-0) = dconjg( dvecdp(itheta,2,0) )
!    dvecdp(itheta,3,0) = - 0.d0 * 0.d0 *inv_sintheta*plm(itheta,1,0) * coef(0) * expimp
     dvecdp(itheta,3,0) = dcmplx( 0.d0,0.d0 )
!    dvecdp(itheta,3,-0) = dconjg( dvecdp(itheta,3,0) )

!! DK DK if m is odd then change the sign of the result
!    if ( mod(m,2)==1 ) then
!       dvec(itheta,1,-m) = - dvec(itheta,1,-m)
!       dvec(itheta,2,-m) = - dvec(itheta,2,-m)
!       dvec(itheta,3,-m) = - dvec(itheta,3,-m)
!       dvecdt(itheta,1,-m) = - dvecdt(itheta,1,-m)
!       dvecdt(itheta,2,-m) = - dvecdt(itheta,2,-m)
!       dvecdt(itheta,3,-m) = - dvecdt(itheta,3,-m)
!       dvecdp(itheta,1,-m) = - dvecdp(itheta,1,-m)
!       dvecdp(itheta,2,-m) = - dvecdp(itheta,2,-m)
!       dvecdp(itheta,3,-m) = - dvecdp(itheta,3,-m)
!    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 1 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (l >= 1) then

     expimp = zexp( dcmplx( 0.d0, phival ) )
     plmdt = x * inv_sintheta * plm(itheta,1,1) + plm(itheta,1,2)
     dvec(itheta,1,1)  = coef(1) * plm(itheta,1,1) * expimp
     dvec(itheta,1,-1) = - dconjg( dvec(itheta,1,1) )
     dvec(itheta,2,1) = coef(1) * plmdt * expimp
     dvec(itheta,2,-1) = - dconjg( dvec(itheta,2,1) )
     dvec(itheta,3,1)  = dcmplx( 0.d0, 1.d0)*inv_sintheta*coef(1)*plm(itheta,1,1)*expimp
     dvec(itheta,3,-1) = - dconjg( dvec(itheta,3,1))

     ! calculate derivatives
     dvecdt(itheta,1,1) = plmdt * coef(1) * expimp
     dvecdt(itheta,1,-1) = - dconjg( dvecdt(itheta,1,1) )
     dvecdt(itheta,2,1) =(-x*inv_sintheta*plmdt+1.d0/(1-x*x)*plm(itheta,1,1)-xl2*plm(itheta,1,1))*coef(1)*expimp
     dvecdt(itheta,2,-1) = - dconjg( dvecdt(itheta,2,1) )
     dvecdt(itheta,3,1) = dcmplx( 0.d0, 1.d0 ) * ( - x / ( 1- x * x ) *plm(itheta,1,1) + inv_sintheta*plmdt)*coef(1)*expimp
     dvecdt(itheta,3,-1) = - dconjg( dvecdt(itheta,3,1) )
     dvecdp(itheta,1,1) = dcmplx( 0.d0, 1.d0 ) * plm(itheta,1,1) * coef(1) * expimp
     dvecdp(itheta,1,-1) = - dconjg( dvecdp(itheta,1,1) )
     dvecdp(itheta,2,1) = dcmplx( 0.d0,1.d0 ) * plmdt * coef(1) * expimp
     dvecdp(itheta,2,-1) = - dconjg( dvecdp(itheta,2,1) )
     dvecdp(itheta,3,1) = - inv_sintheta*plm(itheta,1,1) * coef(1) * expimp
     dvecdp(itheta,3,-1) = - dconjg( dvecdp(itheta,3,1) )

!! DK DK if m is odd then change the sign of the result
!! DK DK moved these minus signs to the expressions in the paragraph above, to reduce memory accesses
!    if ( mod(m,2)==1 ) then
!       dvec(itheta,1,-1) = - dvec(itheta,1,-1)
!       dvec(itheta,2,-1) = - dvec(itheta,2,-1)
!       dvec(itheta,3,-1) = - dvec(itheta,3,-1)
!       dvecdt(itheta,1,-1) = - dvecdt(itheta,1,-1)
!       dvecdt(itheta,2,-1) = - dvecdt(itheta,2,-1)
!       dvecdt(itheta,3,-1) = - dvecdt(itheta,3,-1)
!       dvecdp(itheta,1,-1) = - dvecdp(itheta,1,-1)
!       dvecdp(itheta,2,-1) = - dvecdp(itheta,2,-1)
!       dvecdp(itheta,3,-1) = - dvecdp(itheta,3,-1)
!    endif

  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 2 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (l >= 2) then

     expimp = zexp( dcmplx( 0.d0, 2.d0*phival ) )
     plmdt = 2.d0 * x * inv_sintheta * plm(itheta,1,2) + plm(itheta,1,3)
     dvec(itheta,1,2)  = coef(2) * plm(itheta,1,2) * expimp
     dvec(itheta,1,-2) = dconjg( dvec(itheta,1,2) )
     dvec(itheta,2,2) = coef(2) * plmdt * expimp
     dvec(itheta,2,-2) = dconjg( dvec(itheta,2,2) )
     dvec(itheta,3,2)  = dcmplx( 0.d0, 2.d0)*inv_sintheta*coef(2)*plm(itheta,1,2)*expimp
     dvec(itheta,3,-2) = dconjg( dvec(itheta,3,2))

     ! calculate derivatives
     dvecdt(itheta,1,2) = plmdt * coef(2) * expimp
     dvecdt(itheta,1,-2) = dconjg( dvecdt(itheta,1,2) )
     dvecdt(itheta,2,2) =(-x*inv_sintheta*plmdt+4.d0/(1-x*x)*plm(itheta,1,2)-xl2*plm(itheta,1,2))*coef(2)*expimp
     dvecdt(itheta,2,-2) = dconjg( dvecdt(itheta,2,2) )
     dvecdt(itheta,3,2) = dcmplx( 0.d0, 2.d0 ) * ( - x / ( 1- x * x ) *plm(itheta,1,2) + inv_sintheta*plmdt)*coef(2)*expimp
     dvecdt(itheta,3,-2) = dconjg( dvecdt(itheta,3,2) )
     dvecdp(itheta,1,2) = dcmplx( 0.d0, 2.d0 ) * plm(itheta,1,2) * coef(2) * expimp
     dvecdp(itheta,1,-2) = dconjg( dvecdp(itheta,1,2) )
     dvecdp(itheta,2,2) = dcmplx( 0.d0,2.d0 ) * plmdt * coef(2) * expimp
     dvecdp(itheta,2,-2) = dconjg( dvecdp(itheta,2,2) )
     dvecdp(itheta,3,2) = - 4.d0 *inv_sintheta*plm(itheta,1,2) * coef(2) * expimp
     dvecdp(itheta,3,-2) = dconjg( dvecdp(itheta,3,2) )

!! DK DK if m is odd then change the sign of the result
!    if ( mod(m,2)==1 ) then
!       dvec(itheta,1,-m) = - dvec(itheta,1,-m)
!       dvec(itheta,2,-m) = - dvec(itheta,2,-m)
!       dvec(itheta,3,-m) = - dvec(itheta,3,-m)
!       dvecdt(itheta,1,-m) = - dvecdt(itheta,1,-m)
!       dvecdt(itheta,2,-m) = - dvecdt(itheta,2,-m)
!       dvecdt(itheta,3,-m) = - dvecdt(itheta,3,-m)
!       dvecdp(itheta,1,-m) = - dvecdp(itheta,1,-m)
!       dvecdp(itheta,2,-m) = - dvecdp(itheta,2,-m)
!       dvecdp(itheta,3,-m) = - dvecdp(itheta,3,-m)
!    endif

  endif

! enddo

!! DK DK added this loop (it was previously in the calling program)
  enddo

end subroutine caldvec_for_l_less_than_5_no_store_l

!---------------------------------------------------------------

!! DK DK when l is greater than 5 we can get rid of all the "if" tests in this subroutine to make it much faster
subroutine caldvec_for_l_more_than_5_store_l(l,theta,phi,plm,dvec,dvecdt,dvecdp,theta_n,maxlmax)

  implicit none

  real(kind(0d0)),parameter :: pi=3.1415926535897932d0
  real(kind(0d0)),parameter :: convert_to_radians = pi/180.d0

  integer :: l,m,i,theta_n,itheta,maxlmax

  real(kind(0d0)):: theta(theta_n),phi(theta_n),plm(theta_n,1:3,0:3)
  real(kind(0d0)):: x,fact,thetaval,phival,inv_sintheta
  complex(kind(0d0)) :: dvec(theta_n,1:3,-2:2,0:maxlmax),expimp
  complex(kind(0d0)) :: dvecdt(theta_n,1:3,-2:2,0:maxlmax),dvecdp(theta_n,1:3,-2:2,0:maxlmax)
  real(kind(0d0)):: plmdt,xl2

!! DK DK modified this to precompute the coefficients
  real(kind(0d0)):: coef(0:2)

!! DK DK added this to precompute the coefficients
  do m=0,2
     fact = 1.d0
     if ( m /= 0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef(m) = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
  enddo

!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
  do itheta = 1,theta_n

!! DK DK added this
  thetaval = theta(itheta) * convert_to_radians
  phival = phi(itheta) * convert_to_radians

!! DK DK added this
  inv_sintheta = 1.d0 / dsin(thetaval)

  x = dcos( thetaval )
  xl2 = dble(l) * dble(l+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 0 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     plm(itheta,3,0) = plm(itheta,2,0)
     plm(itheta,2,0) = plm(itheta,1,0)
     plm(itheta,1,0) = (x*dble(2*l-1) * plm(itheta,2,0)-dble(l-1) * plm(itheta,3,0) )/dble(l)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 1 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     plm(itheta,3,1) = plm(itheta,2,1)
     plm(itheta,2,1) = plm(itheta,1,1)
     plm(itheta,1,1) = (x*dble(2*l-1) * plm(itheta,2,1)-dble(l) * plm(itheta,3,1) )/dble(l-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 2 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     plm(itheta,3,2) = plm(itheta,2,2)
     plm(itheta,2,2) = plm(itheta,1,2)
     plm(itheta,1,2) = (x*dble(2*l-1) * plm(itheta,2,2)-dble(l+1) * plm(itheta,3,2) )/dble(l-2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 3 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     plm(itheta,3,3) = plm(itheta,2,3)
     plm(itheta,2,3) = plm(itheta,1,3)
     plm(itheta,1,3) = (x*dble(2*l-1) * plm(itheta,2,3)-dble(l+2) * plm(itheta,3,3) )/dble(l-3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 0 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     expimp = zexp( dcmplx( 0.d0, 0.d0 ) )
     plmdt = plm(itheta,1,1)
     dvec(itheta,1,0,l)  = coef(0) * plm(itheta,1,0) * expimp
     dvec(itheta,2,0,l) = coef(0) * plmdt * expimp
     dvec(itheta,3,0,l)  = dcmplx( 0.d0, 0.d0)

     ! calculate derivatives
     dvecdt(itheta,1,0,l) = plmdt * coef(0) * expimp
     dvecdt(itheta,2,0,l) =(-x*inv_sintheta*plmdt - xl2*plm(itheta,1,0)) *coef(0)*expimp
     dvecdt(itheta,3,0,l) = dcmplx( 0.d0, 0.d0 )
     dvecdp(itheta,1,0,l) = dcmplx( 0.d0, 0.d0 )
     dvecdp(itheta,2,0,l) = dcmplx( 0.d0,0.d0 )
     dvecdp(itheta,3,0,l) = dcmplx( 0.d0,0.d0 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 1 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     expimp = zexp( dcmplx( 0.d0, phival ) )
     plmdt = x * inv_sintheta * plm(itheta,1,1) + plm(itheta,1,2)
     dvec(itheta,1,1,l)  = coef(1) * plm(itheta,1,1) * expimp
     dvec(itheta,1,-1,l) = - dconjg( dvec(itheta,1,1,l) )
     dvec(itheta,2,1,l) = coef(1) * plmdt * expimp
     dvec(itheta,2,-1,l) = - dconjg( dvec(itheta,2,1,l) )
     dvec(itheta,3,1,l)  = dcmplx( 0.d0, 1.d0)*inv_sintheta*coef(1)*plm(itheta,1,1)*expimp
     dvec(itheta,3,-1,l) = - dconjg( dvec(itheta,3,1,l))

     ! calculate derivatives
     dvecdt(itheta,1,1,l) = plmdt * coef(1) * expimp
     dvecdt(itheta,1,-1,l) = - dconjg( dvecdt(itheta,1,1,l) )
     dvecdt(itheta,2,1,l) =(-x*inv_sintheta*plmdt+1.d0/(1-x*x)*plm(itheta,1,1)-xl2*plm(itheta,1,1))*coef(1)*expimp
     dvecdt(itheta,2,-1,l) = - dconjg( dvecdt(itheta,2,1,l) )
     dvecdt(itheta,3,1,l) = dcmplx( 0.d0, 1.d0 ) * ( - x / ( 1- x * x ) *plm(itheta,1,1) + inv_sintheta*plmdt)*coef(1)*expimp
     dvecdt(itheta,3,-1,l) = - dconjg( dvecdt(itheta,3,1,l) )
     dvecdp(itheta,1,1,l) = dcmplx( 0.d0, 1.d0 ) * plm(itheta,1,1) * coef(1) * expimp
     dvecdp(itheta,1,-1,l) = - dconjg( dvecdp(itheta,1,1,l) )
     dvecdp(itheta,2,1,l) = dcmplx( 0.d0,1.d0 ) * plmdt * coef(1) * expimp
     dvecdp(itheta,2,-1,l) = - dconjg( dvecdp(itheta,2,1,l) )
     dvecdp(itheta,3,1,l) = - inv_sintheta*plm(itheta,1,1) * coef(1) * expimp
     dvecdp(itheta,3,-1,l) = - dconjg( dvecdp(itheta,3,1,l) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 2 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     expimp = zexp( dcmplx( 0.d0, 2.d0*phival ) )
     plmdt = 2.d0 * x * inv_sintheta * plm(itheta,1,2) + plm(itheta,1,3)
     dvec(itheta,1,2,l)  = coef(2) * plm(itheta,1,2) * expimp
     dvec(itheta,1,-2,l) = dconjg( dvec(itheta,1,2,l) )
     dvec(itheta,2,2,l) = coef(2) * plmdt * expimp
     dvec(itheta,2,-2,l) = dconjg( dvec(itheta,2,2,l) )
     dvec(itheta,3,2,l)  = dcmplx( 0.d0, 2.d0)*inv_sintheta*coef(2)*plm(itheta,1,2)*expimp
     dvec(itheta,3,-2,l) = dconjg( dvec(itheta,3,2,l))

     ! calculate derivatives
     dvecdt(itheta,1,2,l) = plmdt * coef(2) * expimp
     dvecdt(itheta,1,-2,l) = dconjg( dvecdt(itheta,1,2,l) )
     dvecdt(itheta,2,2,l) =(-x*inv_sintheta*plmdt+4.d0/(1-x*x)*plm(itheta,1,2)-xl2*plm(itheta,1,2))*coef(2)*expimp
     dvecdt(itheta,2,-2,l) = dconjg( dvecdt(itheta,2,2,l) )
     dvecdt(itheta,3,2,l) = dcmplx( 0.d0, 2.d0 ) * ( - x / ( 1- x * x ) *plm(itheta,1,2) + inv_sintheta*plmdt)*coef(2)*expimp
     dvecdt(itheta,3,-2,l) = dconjg( dvecdt(itheta,3,2,l) )
     dvecdp(itheta,1,2,l) = dcmplx( 0.d0, 2.d0 ) * plm(itheta,1,2) * coef(2) * expimp
     dvecdp(itheta,1,-2,l) = dconjg( dvecdp(itheta,1,2,l) )
     dvecdp(itheta,2,2,l) = dcmplx( 0.d0,2.d0 ) * plmdt * coef(2) * expimp
     dvecdp(itheta,2,-2,l) = dconjg( dvecdp(itheta,2,2,l) )
     dvecdp(itheta,3,2,l) = - 4.d0 *inv_sintheta*plm(itheta,1,2) * coef(2) * expimp
     dvecdp(itheta,3,-2,l) = dconjg( dvecdp(itheta,3,2,l) )

  enddo

end subroutine caldvec_for_l_more_than_5_store_l

!---------------------------------------------------------------

subroutine caldvec_for_l_less_than_5_store_l(l,theta,phi,plm,dvec,dvecdt,dvecdp,theta_n,maxlmax)

  implicit none

  real(kind(0d0)),parameter :: pi=3.1415926535897932d0
  real(kind(0d0)),parameter :: convert_to_radians = pi/180.d0

  integer :: l,m,i,theta_n,itheta,maxlmax

  real(kind(0d0)):: theta(theta_n),phi(theta_n),plm(theta_n,1:3,0:3)
  real(kind(0d0)):: x,fact,thetaval,phival,inv_sintheta
  complex(kind(0d0)) :: dvec(theta_n,1:3,-2:2,0:maxlmax),expimp
  complex(kind(0d0)) :: dvecdt(theta_n,1:3,-2:2,0:maxlmax),dvecdp(theta_n,1:3,-2:2,0:maxlmax)
  real(kind(0d0)):: plmdt,xl2

!! DK DK modified this to precompute the coefficients
  real(kind(0d0)):: coef(0:2)

!! DK DK added this to precompute the coefficients
  do m=0,min(l,2)
     fact = 1.d0
     if ( m /= 0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef(m) = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
  enddo

!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
  do itheta = 1,theta_n

!! DK DK added this
  thetaval = theta(itheta) * convert_to_radians
  phival = phi(itheta) * convert_to_radians

!! DK DK added this
  inv_sintheta = 1.d0 / dsin(thetaval)

  x = dcos( thetaval )
  xl2 = dble(l) * dble(l+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 0 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ( l == 0 ) then
     plm(itheta,3,0) = 0.d0
     plm(itheta,2,0) = 0.d0
     plm(itheta,1,0) = 1.d0
  else
     plm(itheta,3,0) = plm(itheta,2,0)
     plm(itheta,2,0) = plm(itheta,1,0)
     if ( l == 1 ) then
        plm(itheta,1,0) = x * plm(itheta,2,0)
     else
        plm(itheta,1,0) = (x*dble(2*l-1) * plm(itheta,2,0)-dble(l-1) * plm(itheta,3,0) )/dble(l)
     endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 1 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (l >= 1) then

  if ( l == 1 ) then
     plm(itheta,3,1) = 0.d0
     plm(itheta,2,1) = 0.d0
     plm(itheta,1,1) = - dsqrt( (1.d0-x)*(1.d0+x) )   !! DK DK computed this analytically instead of with a loop
  else
     plm(itheta,3,1) = plm(itheta,2,1)
     plm(itheta,2,1) = plm(itheta,1,1)
     if ( l == 2 ) then
        plm(itheta,1,1) = x * 3.d0 * plm(itheta,2,1)
     else
        plm(itheta,1,1) = (x*dble(2*l-1) * plm(itheta,2,1)-dble(l) * plm(itheta,3,1) )/dble(l-1)
     endif
  endif

  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 2 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (l >= 2) then

  if ( l == 2 ) then
     plm(itheta,3,2) = 0.d0
     plm(itheta,2,2) = 0.d0
     plm(itheta,1,2) = 3.d0 * dsqrt( (1.d0-x)*(1.d0+x) )**2   !! DK DK computed this analytically instead of with a loop
  else
     plm(itheta,3,2) = plm(itheta,2,2)
     plm(itheta,2,2) = plm(itheta,1,2)
     if ( l == 3 ) then
        plm(itheta,1,2) = x * 5.d0 * plm(itheta,2,2)
     else
        plm(itheta,1,2) = (x*dble(2*l-1) * plm(itheta,2,2)-dble(l+1) * plm(itheta,3,2) )/dble(l-2)
     endif
  endif

  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 3 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (l >= 3) then

  if ( l == 3 ) then
     plm(itheta,3,3) = 0.d0
     plm(itheta,2,3) = 0.d0
     plm(itheta,1,3) = - 15.d0 * dsqrt( (1.d0-x)*(1.d0+x) )**3   !! DK DK computed this analytically instead of with a loop
  else
     plm(itheta,3,3) = plm(itheta,2,3)
     plm(itheta,2,3) = plm(itheta,1,3)
     if ( l == 4 ) then
        plm(itheta,1,3) = x * 7.d0 * plm(itheta,2,3)
     else
        plm(itheta,1,3) = (x*dble(2*l-1) * plm(itheta,2,3)-dble(l+2) * plm(itheta,3,3) )/dble(l-3)
     endif
  endif

  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 0 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     expimp = zexp( dcmplx( 0.d0, 0.d0 ) )
     plmdt = plm(itheta,1,1)
     dvec(itheta,1,0,l)  = coef(0) * plm(itheta,1,0) * expimp
     dvec(itheta,2,0,l) = coef(0) * plmdt * expimp
     dvec(itheta,3,0,l)  = dcmplx( 0.d0, 0.d0)

     ! calculate derivatives
     dvecdt(itheta,1,0,l) = plmdt * coef(0) * expimp
     dvecdt(itheta,2,0,l) =(-x*inv_sintheta*plmdt - xl2*plm(itheta,1,0)) *coef(0)*expimp
     dvecdt(itheta,3,0,l) = dcmplx( 0.d0, 0.d0 )
     dvecdp(itheta,1,0,l) = dcmplx( 0.d0, 0.d0 )
     dvecdp(itheta,2,0,l) = dcmplx( 0.d0,0.d0 )
     dvecdp(itheta,3,0,l) = dcmplx( 0.d0,0.d0 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 1 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (l >= 1) then

     expimp = zexp( dcmplx( 0.d0, phival ) )
     plmdt = x * inv_sintheta * plm(itheta,1,1) + plm(itheta,1,2)
     dvec(itheta,1,1,l)  = coef(1) * plm(itheta,1,1) * expimp
     dvec(itheta,1,-1,l) = - dconjg( dvec(itheta,1,1,l) )
     dvec(itheta,2,1,l) = coef(1) * plmdt * expimp
     dvec(itheta,2,-1,l) = - dconjg( dvec(itheta,2,1,l) )
     dvec(itheta,3,1,l)  = dcmplx( 0.d0, 1.d0)*inv_sintheta*coef(1)*plm(itheta,1,1)*expimp
     dvec(itheta,3,-1,l) = - dconjg( dvec(itheta,3,1,l))

     ! calculate derivatives
     dvecdt(itheta,1,1,l) = plmdt * coef(1) * expimp
     dvecdt(itheta,1,-1,l) = - dconjg( dvecdt(itheta,1,1,l) )
     dvecdt(itheta,2,1,l) =(-x*inv_sintheta*plmdt+1.d0/(1-x*x)*plm(itheta,1,1)-xl2*plm(itheta,1,1))*coef(1)*expimp
     dvecdt(itheta,2,-1,l) = - dconjg( dvecdt(itheta,2,1,l) )
     dvecdt(itheta,3,1,l) = dcmplx( 0.d0, 1.d0 ) * ( - x / ( 1- x * x ) *plm(itheta,1,1) + inv_sintheta*plmdt)*coef(1)*expimp
     dvecdt(itheta,3,-1,l) = - dconjg( dvecdt(itheta,3,1,l) )
     dvecdp(itheta,1,1,l) = dcmplx( 0.d0, 1.d0 ) * plm(itheta,1,1) * coef(1) * expimp
     dvecdp(itheta,1,-1,l) = - dconjg( dvecdp(itheta,1,1,l) )
     dvecdp(itheta,2,1,l) = dcmplx( 0.d0,1.d0 ) * plmdt * coef(1) * expimp
     dvecdp(itheta,2,-1,l) = - dconjg( dvecdp(itheta,2,1,l) )
     dvecdp(itheta,3,1,l) = - inv_sintheta*plm(itheta,1,1) * coef(1) * expimp
     dvecdp(itheta,3,-1,l) = - dconjg( dvecdp(itheta,3,1,l) )

  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! copy for m = 2 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (l >= 2) then

     expimp = zexp( dcmplx( 0.d0, 2.d0*phival ) )
     plmdt = 2.d0 * x * inv_sintheta * plm(itheta,1,2) + plm(itheta,1,3)
     dvec(itheta,1,2,l)  = coef(2) * plm(itheta,1,2) * expimp
     dvec(itheta,1,-2,l) = dconjg( dvec(itheta,1,2,l) )
     dvec(itheta,2,2,l) = coef(2) * plmdt * expimp
     dvec(itheta,2,-2,l) = dconjg( dvec(itheta,2,2,l) )
     dvec(itheta,3,2,l)  = dcmplx( 0.d0, 2.d0)*inv_sintheta*coef(2)*plm(itheta,1,2)*expimp
     dvec(itheta,3,-2,l) = dconjg( dvec(itheta,3,2,l))

     ! calculate derivatives
     dvecdt(itheta,1,2,l) = plmdt * coef(2) * expimp
     dvecdt(itheta,1,-2,l) = dconjg( dvecdt(itheta,1,2,l) )
     dvecdt(itheta,2,2,l) =(-x*inv_sintheta*plmdt+4.d0/(1-x*x)*plm(itheta,1,2)-xl2*plm(itheta,1,2))*coef(2)*expimp
     dvecdt(itheta,2,-2,l) = dconjg( dvecdt(itheta,2,2,l) )
     dvecdt(itheta,3,2,l) = dcmplx( 0.d0, 2.d0 ) * ( - x / ( 1- x * x ) *plm(itheta,1,2) + inv_sintheta*plmdt)*coef(2)*expimp
     dvecdt(itheta,3,-2,l) = dconjg( dvecdt(itheta,3,2,l) )
     dvecdp(itheta,1,2,l) = dcmplx( 0.d0, 2.d0 ) * plm(itheta,1,2) * coef(2) * expimp
     dvecdp(itheta,1,-2,l) = dconjg( dvecdp(itheta,1,2,l) )
     dvecdp(itheta,2,2,l) = dcmplx( 0.d0,2.d0 ) * plmdt * coef(2) * expimp
     dvecdp(itheta,2,-2,l) = dconjg( dvecdp(itheta,2,2,l) )
     dvecdp(itheta,3,2,l) = - 4.d0 *inv_sintheta*plm(itheta,1,2) * coef(2) * expimp
     dvecdp(itheta,3,-2,l) = dconjg( dvecdp(itheta,3,2,l) )

  endif

  enddo

end subroutine caldvec_for_l_less_than_5_store_l

!---------------------------------------------------------------

subroutine now_unused_caldvec_already_plm(theta,phi,plm,dvec,dvecdt,dvecdp,maxlmax,theta_n)

!! DK DK moved the two loops inside the subroutine in order to avoid array copies and to be able to vectorize it

  implicit none

  real(kind(0d0)),parameter :: pi=3.1415926535897932d0
  real(kind(0d0)),parameter :: convert_to_radians = pi/180.d0

  integer, intent(in) :: maxlmax,theta_n
  integer :: l,m,i,itheta
!***************************************************************************************
!***************************************************************************************
!! DK DK l'indice 1:theta_n n'est pas a la bonne place pour vectoriser ici
!! DK DK but it does not matter any more because this subroutine is now unused
!***************************************************************************************
!***************************************************************************************
  real(kind(0d0)), intent(in) :: theta(1:theta_n),phi(1:theta_n),plm(1:3,0:3,1:theta_n,0:maxlmax)
  real(kind(0d0)) :: x,fact,coef,thetaval,phival
  complex(kind(0d0)) :: expimp
  complex(kind(0d0)), intent(out) :: dvec(1:theta_n,1:3,-2:2,0:maxlmax),dvecdt(1:theta_n,1:3,-2:2,0:maxlmax), &
                                        dvecdp(1:theta_n,1:3,-2:2,0:maxlmax)
  real(kind(0d0)) :: plmdt,xl2

  do l = 0,maxlmax
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
     do itheta = 1,theta_n

!! DK DK added this
  thetaval = theta(itheta) * convert_to_radians
  phival = phi(itheta) * convert_to_radians

  x = dcos( thetaval )
  xl2 = dble(l) * dble(l+1)

  do m=0,min0(l,2)
     fact = 1.d0
     if ( m /= 0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     expimp = zexp( dcmplx( 0.d0, dble(m)*phival ) )
     plmdt = dble(m) * x / dsin( thetaval ) * plm(1,m,itheta,l) + plm(1,m+1,itheta,l)

     dvec(itheta,1,m,l)  = coef * plm(1,m,itheta,l) * expimp
     dvec(itheta,1,-m,l) = dconjg( dvec(itheta,1,m,l) )
     dvec(itheta,2,m,l) = coef * plmdt * expimp
     dvec(itheta,2,-m,l) = dconjg( dvec(itheta,2,m,l) )
     dvec(itheta,3,m,l)  = dcmplx( 0.d0, dble(m))/dsin( thetaval )*coef*plm(1,m,itheta,l)*expimp
     dvec(itheta,3,-m,l) = dconjg( dvec(itheta,3,m,l))

     ! calculate derivatives
     dvecdt(itheta,1,m,l) = plmdt * coef * expimp
     dvecdt(itheta,1,-m,l) = dconjg( dvecdt(itheta,1,m,l) )
     dvecdt(itheta,2,m,l) =(-x/dsin(thetaval) * plmdt+dble(m)*dble(m)/(1-x*x)*plm(1,m,itheta,l)-xl2*plm(1,m,itheta,l))*coef*expimp
     dvecdt(itheta,2,-m,l) = dconjg( dvecdt(itheta,2,m,l) )
     dvecdt(itheta,3,m,l) = dcmplx( 0.d0, dble(m) ) * ( - x / ( 1- x * x ) * plm(1,m,itheta,l) + &
                                     1.d0/dsin(thetaval)*plmdt)*coef*expimp
     dvecdt(itheta,3,-m,l) = dconjg( dvecdt(itheta,3,m,l) )
     dvecdp(itheta,1,m,l) = dcmplx( 0.d0, dble(m) ) * plm(1,m,itheta,l) * coef * expimp
     dvecdp(itheta,1,-m,l) = dconjg( dvecdp(itheta,1,m,l) )
     dvecdp(itheta,2,m,l) = dcmplx( 0.d0,dble(m) ) * plmdt * coef * expimp
     dvecdp(itheta,2,-m,l) = dconjg( dvecdp(itheta,2,m,l) )
     dvecdp(itheta,3,m,l) = - dble(m) * dble(m) / dsin(thetaval)*plm(1,m,itheta,l) * coef * expimp
     dvecdp(itheta,3,-m,l) = dconjg( dvecdp(itheta,3,m,l) )

     if ( mod(m,2) == 1 ) then
        dvec(itheta,1,-m,l) = - dvec(itheta,1,-m,l)
        dvec(itheta,2,-m,l) = - dvec(itheta,2,-m,l)
        dvec(itheta,3,-m,l) = - dvec(itheta,3,-m,l)
        dvecdt(itheta,1,-m,l) = - dvecdt(itheta,1,-m,l)
        dvecdt(itheta,2,-m,l) = - dvecdt(itheta,2,-m,l)
        dvecdt(itheta,3,-m,l) = - dvecdt(itheta,3,-m,l)
        dvecdp(itheta,1,-m,l) = - dvecdp(itheta,1,-m,l)
        dvecdp(itheta,2,-m,l) = - dvecdp(itheta,2,-m,l)
        dvecdp(itheta,3,-m,l) = - dvecdp(itheta,3,-m,l)
     endif
  enddo

  enddo
  enddo

end subroutine now_unused_caldvec_already_plm

!---------------------------------------------------------------

subroutine now_unused_clPLM(plm,lmax,theta,theta_n,myrank)

  implicit none

  integer :: itheta, l, theta_n, lmax,m
  real(kind(0d0)) :: plm(1:3,0:3,1:theta_n,0:lmax),theta(1:theta_n)
  real(kind(0d0)) :: tmpthetainrad
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  real(kind(0d0)) :: x, plmtmp(1:3,0:3)
  integer myrank

  do itheta = 1, theta_n
     tmpthetainrad = theta(itheta)/180.d0*pi
     x = cos(tmpthetainrad)
     plmtmp = 0.d0
     do l = 0, lmax
        do m = 0, min0(l,3)
           call calplm(l,m,x,plmtmp(1:3,m))
        enddo
        plm(1:3,0:3,itheta,l) = plmtmp(1:3,0:3)
     enddo
  enddo

end subroutine now_unused_clPLM

!---------------------------------------------------------------

subroutine calplm( l,m,x,plm )
  implicit none
  integer :: l,m,i
  real(kind(0d0)) :: x,plm(1:3),pmm,somx2,fact

  if ((m < 0) .or. (m > l) .or. (dabs(x) > 1.d0)) stop 'bad arguments in calplm'
  if ( l == m ) then
     pmm = 1.d0
     if ( m > 0 ) then
        somx2 = dsqrt( (1.d0-x)*(1.d0+x) )
        fact = 1.d0
        do i=1,m
           pmm = -pmm * fact * somx2
           fact = fact + 2.d0
        enddo
     endif
     plm(3) = 0.d0
     plm(2) = 0.d0
     plm(1) = pmm
  else
     plm(3) = plm(2)
     plm(2) = plm(1)
     if ( l == m+1 ) then
        plm(1) = x * dble(2*m+1) * plm(2)
     else
        plm(1) = (x*dble(2*l-1) * plm(2)-dble(l+m-1) * plm(3) )/dble(l-m)
     endif
  endif

end subroutine calplm

