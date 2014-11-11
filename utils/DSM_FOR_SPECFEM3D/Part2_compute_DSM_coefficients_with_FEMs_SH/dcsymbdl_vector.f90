
!! From Vadim Monteiller, April 2013:

! Je n'ai jamais réussi à comprendre ce qu'il y a dans le système de ces deux routines.
! En dépit des apparences et des commentaires dans le code, le système n'est pas
! symétrique. Nobuaki m'a dit qu'à l'origine ce système était bien
! symétrique mais qu'à la suite des travaux de Ohminato, le système à été
! modifié en ajoutant des coefficients qu'il a trouvé de manière empirique
! qui le rend non symétrique.
!
! Il y a un moment j'avais envisagé d'utiliser MKL, mais je n'ai pas réussi à
! le faire. Roland avait essayé de tester les routines  dcsymbdl0 et dcsbdlv0
! sur un de ses pb, mais il n'obtenait pas la solution voulue, nous avons alors
! conclu que ces routines avaient été modifiées pour les rendre
! opérationnelles sur leur système particulier et que elles ne résolvent
! plus des systèmes symétriques à bandes  diagonales.
!
! Du coup c'est très dur de savoir ce qu'il y a dans les tableaux, même en
! discutant avec Nobuaki j'ai pas du tout compris comment étaient stockées
! les différentes diagonales, j'ai juste cru comprendre que c'était pas la
! même convention que les routines de MKL.  Du coup, je ne sais pas si ça
! vaut le coup de regarder çà à fond, surtout que l'essentiel du coût
! était dans le relecture des coefs et sur ce point-là on a beaucoup
! gagné. Je pense qu'il vaut mieux qu'on regarde les I/O parallèles.

subroutine dcsymbdl0_m_equals_3_nn_equals_6(a,n,eps,z,w,ier)

  !**********************************************************************
  !                                                                     *
  !  Gauss method for a symmetric band matrix with a short width.       *
  !  this routine is for a vector computer.                             *
  !                                                                     *
  !  parameters:                                                        *
  !   on entry:                                                         *
  !     a      the array with dimension (m+1)*n which contains          *
  !            the lower band matrix in row-wise.                       *
  !     m      the half band width of the matrix a, excluding diagonals.*
  !     n      the order of the matrix a.                               *
  !     nn     the order of working arrays l, li, and lj.               *
  !     eps    the tolerance for pivotal elements.                      *
  !   on return:                                                        *
  !     a      the lower triangular matrix l after decomposition.       *
  !     ier    error code. if ier=0, normal return.                     *
  !   others: working parameters.                                       *
  !                                                                     *
  ! Copyright: Fumiko Nagahori 1 Sep 1991, version 1                    *
  !  Dimitri Komatitsch April 2013, version 2: optimized and vectorized *
  !                                                                     *
  !**********************************************************************

  implicit none

!! DK DK made a specific version in order to be able to fully optimize

  integer :: n
  integer :: ier

  real(kind(0d0)) :: eps
  complex(kind(0d0)) :: a(4*n),z(4),w(4)
  integer :: i,j,k,nk,nkk,nki
  complex(kind(0d0)) :: pivot

  complex(kind(0d0)) :: z2,z3,z4
  complex(kind(0d0)) :: w2,w3,w4

! ij = 0
! do i=1,3
!   ij = ij + 1
!   li(ij) = i + 1
!   lj(ij) = 2
!   l(ij) = 3 * (i-1) + 1
!   do j=i+1,3
!     ij = ij + 1
!     l(ij) = l(ij-1) + 4
!     li(ij) = li(ij-1) + 1
!     lj(ij) = lj(ij-1) + 1
!   enddo
! enddo
!
!! DK DK the above loops give the following values:
! integer, parameter :: l1  =  1, l2  =  5, l3  =  9, l4  =  4, l5  =  8, l6  =  7
! integer, parameter :: li1  = 2, li2  = 3, li3  = 4, li4  = 3, li5  = 4, li6  = 4
! integer, parameter :: lj1  = 2, lj2  = 3, lj3  = 4, lj4  = 2, lj5  = 3, lj6  = 2

  ier = 0

!! DK DK this is the main loop on the whole matrix except its last two elements, which will be handled separately later;
!! DK DK it cannot be vectorized because of an intrinsic vector dependence (since it implements Gaussian elimination)
!IBM* NOVECTOR
!DIR$ NOVECTOR
  do k=1,n-3

!! DK DK
    nk = 4*k

! error if the pivot is almost zero
    if (zabs(a(nk)) < eps) then
      print *, '(subr. symbdl) singular at step = ', k
      ier = 1
      return
    endif

    pivot = (1.0d0,0.0d0) / a(nk)

    z2 = - a(nk+3)
    w2 =   a(nk+3) * pivot

    z3 = - a(nk+6)
    w3 =   a(nk+6) * pivot

    z4 = - a(nk+9)
    w4 =   a(nk+9) * pivot

!! DK DK these are probably the formulas from Ohminato mentioned by Nobuaki
    a(nk+3)  = w2
    a(nk+4)  = a(nk+4)  + w2 * z2
!! DK DK nk+5 is unchanged here for some reason
    a(nk+6)  = w3
    a(nk+7)  = a(nk+7)  + w2 * z3
    a(nk+8)  = a(nk+8)  + w3 * z3
    a(nk+9)  = w4
    a(nk+10) = a(nk+10) + w2 * z4
    a(nk+11) = a(nk+11) + w3 * z4
    a(nk+12) = a(nk+12) + w4 * z4
  enddo

!! DK DK these are the additional loops handling the last two elements;
!! DK DK these nested loops do not need to be optimized because they apply to the last two indices only
  do k=n-2,n-1
     nk = 4 * (k-1) + 1
     nkk = 4 * k - 1

     if (zabs(a(nk+3)) < eps) then
        print *, '(subr. symbdl) singular at step = ', k
        ier = 1
        return
     endif

     pivot = (1.0d0,0.0d0) / a(nk+3)

     do j=2,n-k+1
       z(j) = - a(nk+3*j)
       w(j) =   a(nk+3*j) * pivot
       a(nk+3*j) = w(j)
       nki = nkk + 3 * (j-1)
       do i=2,j
         a(nki+i) = a(nki+i) + w(j) * z(i)
       enddo
     enddo
  enddo

end subroutine dcsymbdl0_m_equals_3_nn_equals_6

!--------------------------------------------------------------

subroutine dcsbdlv0_m_equals_3(a,b,n,z)

  !**********************************************************************
  !                                                                     *
  !  Gauss method for a symmetric band matrix with a short width.       *
  !  this routine is for a vector computer.                             *
  !                                                                     *
  !  parameters:                                                        *
  !   on entry:                                                         *
  !     a      the array with dimension (m+1),n which contains          *
  !            the left hand side lower band matrix.                    *
  !     m      the half band width of the matrix a, excluding diagonals.*
  !     n      the order of the matrix a.                               *
  !     b      right hand side vector                                   *
  !   on return:                                                        *
  !     b      solution.                                                *
  !   others: working parameters.                                       *
  ! Copyright: Fumiko Nagahori 1 Sep 1991, version 1                    *
  !  Dimitri Komatitsch April 2013, version 2: optimized and vectorized *
  !                                                                     *
  !**********************************************************************

  implicit none

!! DK DK made a specific version in order to be able to fully optimize and vectorize

  integer :: n
  complex(kind(0d0)) :: a(4,n),b(n),z(n)
  integer :: j,k,i1,j1
  complex(kind(0d0)) :: sum_values

! forward substitution

    z(1) = b(1)
    z(2) = b(2) - a(2,2) * z(1)

    do j=3,n
      if (j > 4) then
        i1 = 1
      else
        i1 = 5 - j
      endif
      sum_values = (0.0d0,0.0d0)
!! DK DK i1 has a big negative value here in most cases
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
      do k=i1,3
        sum_values = sum_values + a(k,j) * z(j-4+k)
      enddo
      z(j) = b(j) - sum_values
    enddo

!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
    do j=1,n
      z(j) = z(j) / a(4,j)
    enddo

    b(n) = z(n)
    b(n-1) = z(n-1) - a(3,n) * z(n)

    do j=3,n
      j1 = n - j + 1
      i1 = 1
      if (j < 4) i1 = 5 - j
      sum_values = (0.0d0,0.0d0)
!! DK DK i1 has a big negative value here in most cases
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
      do k=i1,3
        sum_values = sum_values + a(k,4-k+j1) * b(4-k+j1)
      enddo
      b(j1) = z(j1) - sum_values
    enddo

end subroutine dcsbdlv0_m_equals_3

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

!! DK DK no need to optimize the version below, which is called only once, for l == 0

subroutine dcsymbdl0_m_equals_1_nn_equals_1(a,n,eps,z,w,l,li,lj,ier)

  !**********************************************************************
  !                                                                     *
  !  Gauss method for a symmetric band matrix with a short width.       *
  !  this routine is for a vector computer.                             *
  !                                                                     *
  !  parameters:                                                        *
  !   on entry:                                                         *
  !     a      the array with dimension (m+1)*n which contains          *
  !            the lower band matrix in row-wise.                       *
  !     m      the half band width of the matrix a, excluding diagonals.*
  !     n      the order of the matrix a.                               *
  !     nn     the order of working arrays l, li, and lj.               *
  !     eps    the tolerance for pivotal elements.                      *
  !   on return:                                                        *
  !     a      the lower triangular matrix l after decomposition.       *
  !     ier    error code. if ier=0, normal return.                     *
  !   others: working parameters.                                       *
  !                                                                     *
  ! Copyright: Fumiko Nagahori 1 Sep 1991, version 1                    *
  !  Dimitri Komatitsch April 2013, version 2: optimized and vectorized *
  !                                                                     *
  !**********************************************************************

  implicit none

!! WY made a specific version in order to be able to fully optimize
  integer, parameter :: m = 1, nn = 1, mm = (m+1) * m / 2

  integer :: n
  integer :: l(nn),li(nn),lj(nn),ier
  real(kind(0d0)) :: eps
  complex(kind(0d0)) :: a((m+1)*n),z(m+1),w(m+1)
  integer :: i,j,k,ij,kk,nk,nkk,nki
  complex(kind(0d0)) :: pivot


! mm = 1 at the moment
  ier = 0
  ij = 0

!  do i=1,m  !m = 1
    ij = ij + 1
    li(ij) = 2
    lj(ij) = 2
    l(ij) = 1 !m = 1
!    do j=i+1,m  !m = 1
!      ij = ij + 1
!      l(ij) = l(ij-1) + m + 1
!      li(ij) = li(ij-1) + 1
!      lj(ij) = lj(ij-1) + 1
!    enddo
!  enddo
!WY because m = 1, so the inner cycle don't run.
!WY and at the moment, the ij = 1. work array li, lj and l are only scalar.


  do k=1,n-1   !m = 1
    nk = 2 * (k-1) + 1

    if (zabs(a(nk+1)) < eps) then
      print *, '(subr. symbdl) singular at step = ', k
      ier = 1
      return
    endif

    pivot = (1.0d0,0.0d0) / a(nk+1)  !m = 1

!    do j=2,m+1  ! m = 1
      j = 2
      z(j) = - a(nk+2)
      w(j) =   a(nk+2) * pivot
      a(nk+2) = w(j)
!    enddo

    kk = nk + 2

    ! vortion vec
!    do i=1,mm !mm = 1
      a(kk+1) = a(kk+1) + w(2) * z(2)
!    enddo

  enddo

!  do k=n-m+1,n-1  ! because m =1 ,so the loop don't run
!     nk = (m+1) * (k-1) + 1
!     nkk = (m+1) * k - 1

!     if (zabs(a(nk+m)) < eps) then
!        print *, '(subr. symbdl) singular at step = ', k
!        ier = 1
!        return
!     endif

!     pivot = (1.0d0,0.0d0) / a(nk+m)

!     do j=2,n-k+1
!       z(j) = - a(nk+m*j)
!       w(j) =   a(nk+m*j) * pivot
!       a(nk+m*j) = w(j)
!       nki = nkk + m * (j-1)
!       do i=2,j
!         a(nki+i) = a(nki+i) + w(j) * z(i)
!       enddo
!     enddo
!  enddo

end subroutine dcsymbdl0_m_equals_1_nn_equals_1

!--------------------------------------------------------------

subroutine dcsbdlv0_m_equals_1(a,b,n,z)

  !**********************************************************************
  !                                                                     *
  !  Gauss method for a symmetric band matrix with a short width.       *
  !  this routine is for a vector computer.                             *
  !                                                                     *
  !  parameters:                                                        *
  !   on entry:                                                         *
  !     a      the array with dimension (m+1),n which contains          *
  !            the left hand side lower band matrix.                    *
  !     m      the half band width of the matrix a, excluding diagonals.*
  !     n      the order of the matrix a.                               *
  !     b      right hand side vector                                   *
  !   on return:                                                        *
  !     b      solution.                                                *
  !   others: working parameters.                                       *
  ! Copyright: Fumiko Nagahori 1 Sep 1991, version 1                    *
  !  Dimitri Komatitsch April 2013, version 2: optimized and vectorized *
  !                                                                     *
  !**********************************************************************

  implicit none

  integer, parameter :: m = 1

  integer :: n
  complex(kind(0d0)) :: a(2,n),b(n),z(n)
  integer :: j

! forward substitution

    z(1) = b(1)

!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
    do j=2,n
      z(j) = b(j) - a(1,j) * z(j-1)
    enddo

!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
    do j=1,n
      z(j) = z(j) / a(2,j)
    enddo
    b(n) = z(n)

!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
    do j=1,n-1
      b(n-j) = z(n-j) - a(1,n-j+1) * b(n-j+1)
    enddo

end subroutine dcsbdlv0_m_equals_1

