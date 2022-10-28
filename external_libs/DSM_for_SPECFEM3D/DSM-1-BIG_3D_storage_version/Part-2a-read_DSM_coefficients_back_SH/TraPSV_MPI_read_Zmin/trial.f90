subroutine calbvecphi0_stock( l,thetadeg,plmDir,bvec,bvecdt,bvecdp)

  implicit none
  character(*) :: plmDir
  character(120) :: coutfile
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0
  integer  :: l,m,i,j
  real(kind(0d0)) :: theta,thetadeg,x,plm(1:3,0:3),fact,coef
  complex(kind(0d0)) :: bvec(1:3,-2:2)
  complex(kind(0d0)) :: bvecdt(1:3,-2:2),bvecdp(1:3,-2:2)
  real(kind(0d0)) :: plmdt,xl2


  write(coutfile, '(I8, ".","PLM")') int(thetadeg*100000.d0)
  do j = 1,8
     if (coutfile(j:j) == ' ')coutfile(j:j) = '0'
  enddo

  coutfile = plmDir//"/"//coutfile

  open(1,file=coutfile,status='old',form='unformatted',access='direct', &
       recl=kind(0d0)*12)
  read(1,rec=l+1)plm(1,0),plm(1,1),plm(1,2),plm(1,3), &
       plm(2,0),plm(2,1),plm(2,2),plm(2,3), &
       plm(3,0),plm(3,1),plm(3,2),plm(3,3)
  close(1)

  theta =  thetadeg/180.d0*pi

  x = dcos( theta )
  xl2 = dble(l) * dble(l+1)

  do m=0,min0(l,2)
     fact = 1.d0
     if ( m /= 0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     plmdt = dble(m) * x / sin( theta ) * plm(1,m) + plm(1,m+1)
     bvec(1,m)  = dcmplx( 0.d0 )
     bvec(1,-m) = dcmplx( 0.d0 )
     bvec(2,m)  = dcmplx( 0.d0, dble(m) ) / dsin( theta) * coef * plm(1,m)
     bvec(2,-m) = dcmplx(conjg( bvec(2,m)) )
     bvec(3,m) = - coef * plmdt
     bvec(3,-m) = dcmplx(conjg( bvec(3,m) ))

     ! calculate derivatives
     bvecdt(1,m)  = dcmplx( 0.d0 )
     bvecdt(1,-m) = dcmplx( 0.d0 )
     bvecdt(2,m)  = dcmplx( 0.d0, dble(m) ) * ( plmdt / dsin(theta) &
          - x / ( 1 - x * x ) * plm(1,m) ) * coef
     bvecdt(2,-m) = dcmplx( conjg( bvecdt(2,m) ))
     bvecdt(3,m) = ( x / dsin(theta) * plmdt - dble(m) * dble(m)/(1-x*x) *plm(1,m) &
          &           + xl2 * plm(1,m) ) * coef
     bvecdt(3,-m) = dcmplx(conjg( bvecdt(3,m)) )
     bvecdp(1,m)  = dcmplx( 0.d0 )
     bvecdp(1,-m) = dcmplx( 0.d0 )
     bvecdp(2,m)  = - dble(m) * dble(m) / dsin(theta) * plm(1,m) * coef
     bvecdp(2,-m) = dcmplx(conjg( bvecdp(2,m)) )
     bvecdp(3,m)  = - dcmplx( 0.d0, dble(m) ) * plmdt * coef
     bvecdp(3,-m) = dcmplx(conjg( bvecdp(3,m)) )

     if ( mod(m,2) == 1 ) then
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
        bvecdt(2,-m) = - bvecdt(2,-m)
        bvecdt(3,-m) = - bvecdt(3,-m)
        bvecdp(2,-m) = - bvecdp(2,-m)
        bvecdp(3,-m) = - bvecdp(3,-m)
     endif
  enddo

  return
end subroutine calbvecphi0_stock



subroutine calbvecphi0( l,theta,plm,bvec,bvecdt,bvecdp)

  implicit none
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0
  integer  :: l,m,i
  real(kind(0d0)) :: theta,x,plm(1:3,0:3),fact,coef
  complex(kind(0d0)) :: bvec(1:3,-2:2)
  complex(kind(0d0)) :: bvecdt(1:3,-2:2),bvecdp(1:3,-2:2)
  real(kind(0d0)) :: plmdt,xl2

  x = dcos( theta )
  xl2 = dble(l) * dble(l+1)
  do m=0,min0(l,3)
     call calplm( l,m,x,plm(1:3,m))
  enddo

  do m=0,min0(l,2)
     fact = 1.d0
     if ( m /= 0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     plmdt = dble(m) * x / sin( theta ) * plm(1,m) + plm(1,m+1)
     bvec(1,m)  = dcmplx( 0.d0 )
     bvec(1,-m) = dcmplx( 0.d0 )
     bvec(2,m)  = dcmplx( 0.d0, dble(m) ) / dsin( theta) * coef * plm(1,m)
     bvec(2,-m) = dcmplx(conjg( bvec(2,m)) )
     bvec(3,m) = - coef * plmdt
     bvec(3,-m) = dcmplx(conjg( bvec(3,m) ))

     ! calculate derivatives
     bvecdt(1,m)  = dcmplx( 0.d0 )
     bvecdt(1,-m) = dcmplx( 0.d0 )
     bvecdt(2,m)  = dcmplx( 0.d0, dble(m) ) * ( plmdt / dsin(theta) &
          - x / ( 1 - x * x ) * plm(1,m) ) * coef
     bvecdt(2,-m) = dcmplx( conjg( bvecdt(2,m) ))
     bvecdt(3,m) = ( x / dsin(theta) * plmdt - dble(m) * dble(m)/(1-x*x) *plm(1,m) &
          &           + xl2 * plm(1,m) ) * coef
     bvecdt(3,-m) = dcmplx(conjg( bvecdt(3,m)) )
     bvecdp(1,m)  = dcmplx( 0.d0 )
     bvecdp(1,-m) = dcmplx( 0.d0 )
     bvecdp(2,m)  = - dble(m) * dble(m) / dsin(theta) * plm(1,m) * coef
     bvecdp(2,-m) = dcmplx(conjg( bvecdp(2,m)) )
     bvecdp(3,m)  = - dcmplx( 0.d0, dble(m) ) * plmdt * coef
     bvecdp(3,-m) = dcmplx(conjg( bvecdp(3,m)) )

     if ( mod(m,2) == 1 ) then
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
        bvecdt(2,-m) = - bvecdt(2,-m)
        bvecdt(3,-m) = - bvecdt(3,-m)
        bvecdp(2,-m) = - bvecdp(2,-m)
        bvecdp(3,-m) = - bvecdp(3,-m)
     endif
  enddo



  return
end subroutine calbvecphi0



subroutine calplm( l,m,x,plm )
  implicit none
  integer :: l,m,i
  real(kind(0d0)) :: x,plm(1:3),pmm,somx2,fact

  if ((m < 0) .or. (m > l) .or. (dabs(x) > 1.d0)) pause 'bad arguments'
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
        !print *, l,m,x
        plm(1) = (x*dble(2*l-1) * plm(2)-dble(l+m-1) * plm(3) )/dble(l-m)
        !print *, plm(1)
     endif
  endif


end subroutine calplm


subroutine calbveczero( l,bvec )

  implicit none
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0

  integer  :: l,m,i
  real(kind(0d0)) :: fact,coef
  complex(kind(0d0)) :: bvec(1:3,-2:2)
  real(kind(0d0)) :: xl2


  xl2 = dble(l) * dble(l+1)

  do m=0,min0(l,1)
     fact = 1.d0
     if ( m /= 0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     bvec(1,m)  = dcmplx( 0.d0 )
     bvec(1,-m) = dcmplx( 0.d0 )
     bvec(2,m)  = dcmplx( 0.d0, dble(m)) *  xl2 *coef / 2.d0
     bvec(2,-m) = dcmplx(conjg( bvec(2,m)) )
     bvec(3,m) =  -dcmplx(dble(m),0.d0) * xl2 * coef / 2.d0
     bvec(3,-m) = dcmplx(conjg( bvec(3,m) ))

     if ( mod(m,2) == 1 ) then
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
     endif
  enddo



  return
end subroutine calbveczero

subroutine calbvec( l,theta,phi,plm,bvec,bvecdt,bvecdp )
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c Evaluating the value of toroidal harmonics (fully normalized)
      !c at each station whose latitude and longitude are theta and phi.
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  real(kind(0d0)),parameter :: pi=3.1415926535897932d0

  integer ::  l,m,i
  real(kind(0d0)):: theta,phi,x,plm(3,0:3),fact,coef
  complex(kind(0d0)) :: bvec(3,-2:2),expimp
  complex(kind(0d0)) :: bvecdt(3,-2:2),bvecdp(3,-2:2)
  real(kind(0d0)) :: plmdt,xl2

 x = dcos( theta )
 xl2 = dble(l) * dble(l+1)
 do m=0,min0(l,3)
   call calplm( l,m,x,plm(1,m) )
 enddo

 do m=0,min0(l,2)
  fact = 1.d0
   if ( m /= 0 ) then
    do i=l-m+1,l+m
     fact = fact * dble(i)
    enddo
   endif

  coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
  expimp = cdexp( dcmplx( 0.d0, dble(m)*phi ) )
  plmdt = dble(m) * x / dsin( theta ) * plm(1,m) + plm(1,m+1)
  bvec(1,m)  = dcmplx( 0.d0 )
  bvec(1,-m) = dcmplx( 0.d0 )
  bvec(2,m)  = dcmplx( 0.d0, dble(m) ) / dsin( theta ) * coef * plm(1,m) * expimp
  bvec(2,-m) = dconjg( bvec(2,m) )
  bvec(3,m) = - coef * plmdt * expimp
  bvec(3,-m) = dconjg( bvec(3,m) )
  ! calculate derivatives
  bvecdt(1,m)  = dcmplx( 0.d0 )
  bvecdt(1,-m) = dcmplx( 0.d0 )
  bvecdt(2,m)  = dcmplx( 0.d0, dble(m) )* ( plmdt / dsin(theta)- x/(1-x*x)*plm(1,m))*coef*expimp
  bvecdt(2,-m) = dconjg( bvecdt(2,m) )
  bvecdt(3,m) = (x/dsin(theta)*plmdt-dble(m)*dble(m)/(1-x*x)*plm(1,m)+xl2*plm(1,m))*coef*expimp
  bvecdt(3,-m) = dconjg( bvecdt(3,m) )
  bvecdp(1,m)  = dcmplx( 0.d0 )
  bvecdp(1,-m) = dcmplx( 0.d0 )
  bvecdp(2,m)  = - dble(m) * dble(m) / dsin(theta) * plm(1,m)*coef * expimp
  bvecdp(2,-m) = dconjg( bvecdp(2,m) )
  bvecdp(3,m)  = - dcmplx( 0.d0, dble(m) ) * plmdt * coef * expimp
  bvecdp(3,-m) = dconjg( bvecdp(3,m) )

  if ( mod(m,2) == 1 ) then
   bvec(2,-m) = - bvec(2,-m)
   bvec(3,-m) = - bvec(3,-m)
   bvecdt(2,-m) = - bvecdt(2,-m)
   bvecdt(3,-m) = - bvecdt(3,-m)
   bvecdp(2,-m) = - bvecdp(2,-m)
   bvecdp(3,-m) = - bvecdp(3,-m)
  endif
 enddo
 return
end subroutine calbvec


subroutine calbvec_vector( l,theta,phi,plm,bvec,bvecdt,bvecdp,theta_n)
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c Evaluating the value of toroidal harmonics (fully normalized)
  !c at each station whose latitude and longitude are theta and phi.
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  real(kind(0d0)),parameter :: pi=3.1415926535897932d0

  integer ::  l,m,i,theta_n,itheta
  real(kind(0d0)):: theta(1:theta_n),phi(1:theta_n),x,plm(1:theta_n,1:3,0:3),fact,coef(0:2),inv_sintheta
  complex(kind(0d0)) :: bvec(1:theta_n,1:3,-2:2),expimp
  complex(kind(0d0)) :: bvecdt(1:theta_n,1:3,-2:2),bvecdp(theta_n,1:3,-2:2)
  real(kind(0d0)) :: plmdt,xl2

!write(*,*) 'l=',l,'theta_n is',theta_n

  xl2 = dble(l) * dble(l+1)

 if (l < 5 ) then
   do itheta = 1,theta_n
   x = dcos( theta(itheta) )

    do m=0,min0(l,3)
       call calplm( l,m,x,plm(itheta,1:3,m) )
    enddo
   enddo
 else
   do itheta = 1,theta_n
   x = dcos( theta(itheta) )
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

   enddo
 endif


! plm computation is over.

  do m=0,min0(l,2)
     fact = 1.d0
     if ( m /= 0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef(m) = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
  enddo
! unify cofficients computation.


 do itheta = 1,theta_n
 x = dcos( theta(itheta) )  ! please noticing don't miss this x expression.
 inv_sintheta = 1.d0 / dsin(theta(itheta)) !pre-computing the inverse of sin(theta)
! this is for m==0, for any l value we need to calculate thses.
     expimp = cdexp( dcmplx( 0.d0,0.d0 ))
     plmdt = plm(itheta,1,1)
     bvec(itheta,1,0)  = dcmplx( 0.d0 ,0.d0)
     bvec(itheta,2,0)  = dcmplx( 0.d0, 0.d0)
     bvec(itheta,3,0) = - coef(0) * plmdt * expimp
     ! calculate derivatives
     bvecdt(itheta,1,0)  = dcmplx( 0.d0, 0.d0)
     bvecdt(itheta,2,0)  = dcmplx( 0.d0, 0.d0)
     bvecdt(itheta,3,0) = (x*inv_sintheta*plmdt + xl2*plm(itheta,1,0))*coef(0)*expimp
     bvecdp(itheta,1,0)  = dcmplx( 0.d0, 0.d0)
     bvecdp(itheta,2,0)  = dcmplx( 0.d0, 0.d0)
     bvecdp(itheta,3,0)  = dcmplx( 0.d0, 0.d0)
! this is for m==1, for any l >= 1
  if (l >= 1) then

     expimp = cdexp( dcmplx( 0.d0, phi(itheta) ) )
     plmdt =  x *inv_sintheta  * plm(itheta,1,1) + plm(itheta,1,2)
     bvec(itheta,1,1)  = dcmplx( 0.d0 )
     bvec(itheta,1,-1) = dcmplx( 0.d0 )
     bvec(itheta,2,1)  = dcmplx( 0.d0, 1.d0 ) *inv_sintheta * coef(1) * plm(itheta,1,1) * expimp
     bvec(itheta,2,-1) = - dconjg( bvec(itheta,2,1) )
     bvec(itheta,3,1) = - coef(1) * plmdt * expimp
     bvec(itheta,3,-1) = - dconjg( bvec(itheta,3,1) )
     ! calculate derivatives
     bvecdt(itheta,1,1)  = dcmplx( 0.d0 )
     bvecdt(itheta,1,-1) = dcmplx( 0.d0 )
     bvecdt(itheta,2,1)  = dcmplx( 0.d0, 1.d0 )* ( plmdt *inv_sintheta - x/(1-x*x)*plm(itheta,1,1))*coef(1)*expimp
     bvecdt(itheta,2,-1) = - dconjg( bvecdt(itheta,2,1) )
     bvecdt(itheta,3,1) = (x*inv_sintheta*plmdt-1.d0/(1-x*x)*plm(itheta,1,1)+xl2*plm(itheta,1,1))*coef(1)*expimp
     bvecdt(itheta,3,-1) = - dconjg( bvecdt(itheta,3,1) )
     bvecdp(itheta,1,1)  = dcmplx( 0.d0 )
     bvecdp(itheta,1,-1) = dcmplx( 0.d0 )
     bvecdp(itheta,2,1)  = - inv_sintheta * plm(itheta,1,1)*coef(1) * expimp
     bvecdp(itheta,2,-1) = - dconjg( bvecdp(itheta,2,1) )
     bvecdp(itheta,3,1)  = - dcmplx( 0.d0, 1.d0 ) * plmdt * coef(1) * expimp
     bvecdp(itheta,3,-1) = - dconjg( bvecdp(itheta,3,1) )

  endif
! this is for m==2, for any l >= 2

  if (l >= 2) then
     expimp = cdexp( dcmplx( 0.d0, 2.d0*phi(itheta) ) )
     plmdt = 2.d0 * x *inv_sintheta * plm(itheta,1,2) + plm(itheta,1,3)
     bvec(itheta,1,2)  = dcmplx( 0.d0 )
     bvec(itheta,1,-2) = dcmplx( 0.d0 )
     bvec(itheta,2,2)  = dcmplx( 0.d0, 2.d0 ) *inv_sintheta * coef(2) * plm(itheta,1,2) * expimp
     bvec(itheta,2,-2) = dconjg( bvec(itheta,2,2) )
     bvec(itheta,3,2) = - coef(2) * plmdt * expimp
     bvec(itheta,3,-2) = dconjg( bvec(itheta,3,2) )
     ! calculate derivatives
     bvecdt(itheta,1,2)  = dcmplx( 0.d0 )
     bvecdt(itheta,1,-2) = dcmplx( 0.d0 )
     bvecdt(itheta,2,2)  = dcmplx( 0.d0, 2.d0 )* ( plmdt*inv_sintheta - x/(1-x*x)*plm(itheta,1,2))*coef(2)*expimp
     bvecdt(itheta,2,-2) = dconjg( bvecdt(itheta,2,2) )
     bvecdt(itheta,3,2) = (x*inv_sintheta*plmdt-4.d0/(1-x*x)*plm(itheta,1,2)+xl2*plm(itheta,1,2))*coef(2)*expimp
     bvecdt(itheta,3,-2) = dconjg( bvecdt(itheta,3,2) )
     bvecdp(itheta,1,2)  = dcmplx( 0.d0 )
     bvecdp(itheta,1,-2) = dcmplx( 0.d0 )
     bvecdp(itheta,2,2)  = - 4.d0 *inv_sintheta * plm(itheta,1,2)*coef(2) * expimp
     bvecdp(itheta,2,-2) = dconjg( bvecdp(itheta,2,2) )
     bvecdp(itheta,3,2)  = - dcmplx( 0.d0, 2.d0 ) * plmdt * coef(2) * expimp
     bvecdp(itheta,3,-2) = dconjg( bvecdp(itheta,3,2) )

  endif


!     if ( mod(m,2)==1 ) then
!        bvec(2,-m) = - bvec(2,-m)
!        bvec(3,-m) = - bvec(3,-m)
!        bvecdt(2,-m) = - bvecdt(2,-m)
!        bvecdt(3,-m) = - bvecdt(3,-m)
!        bvecdp(2,-m) = - bvecdp(2,-m)
!        bvecdp(3,-m) = - bvecdp(3,-m)
!     endif

 enddo
  return
end subroutine calbvec_vector




subroutine calbvec_vector_test( l,theta,phi,plm,bvec,bvecdt,bvecdp,theta_n)
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c Evaluating the value of toroidal harmonics (fully normalized)
  !c at each station whose latitude and longitude are theta and phi.
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  real(kind(0d0)),parameter :: pi=3.1415926535897932d0

  integer ::  l,m,i,theta_n,itheta
  real(kind(0d0)):: theta(theta_n),phi(theta_n),x,plm(theta_n,1:3,0:3),fact,coef(0:2),inv_sintheta
  complex(kind(0d0)) :: bvec(theta_n,1:3,-2:2),expimp
  complex(kind(0d0)) :: bvecdt(theta_n,1:3,-2:2),bvecdp(theta_n,1:3,-2:2)
  real(kind(0d0)) :: plmdt,xl2

!write(*,*) 'l=',l,'theta_n is',theta_n

  xl2 = dble(l) * dble(l+1)
 do itheta = 1,theta_n
  x = dcos( theta(itheta) )

  if (l < 5 ) then
  write(*,*) 'if l < 5',l
   do m=0,min0(l,3)
    write(*,*) 'when l < 5,the m =',m
      call calplm( l,m,x,plm(itheta,1,m) )

    if (itheta == 6)  write(*,*)  'itheta th plm components are',m,plm(itheta,:,m)

   enddo

   if (itheta == 6) then
!   write(*,*) itheta,'th station',itheta,theta(itheta),x,l,m
!   write(*,*) 'itheta th plm components are:1',plm(itheta,1,:)
!   write(*,*) 'itheta th plm components are:2',plm(itheta,2,:)
!   write(*,*) 'itheta th plm components are:3',plm(itheta,3,:)
   write(*,*) '\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
   endif


  else
  write(*,*) 'if l >= 5',l
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

  endif
 enddo

 itheta = 6
 write(*,*) 'itheta th plm components are:0',plm(itheta,:,0)
 write(*,*) 'itheta th plm components are:1',plm(itheta,:,1)
 write(*,*) 'itheta th plm components are:2',plm(itheta,:,2)
 write(*,*) 'itheta th plm components are:3',plm(itheta,:,3)
! plm computation is over.

  do m=0,min0(l,2)
     fact = 1.d0
     if ( m /= 0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef(m) = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
  enddo
! unify cofficients computation.

 write(*,*) 'coef are:',coef

end subroutine calbvec_vector_test
