!
!
!   ROUTINES POUR FAIRE DES ROTATIONS 3D ET DIVERS CHANGEMENTS DE REPERES
!
! Vadim Monteiller Mars 2013
!
!-------------------------------------------------------------------------------
! matrice de rotation 3D d'axe "axe" et d'angle theta (degres)
! cette matrice est en complexe
subroutine rotation_matrix(R,axe,theta)
  implicit none
  double precision axe(3),theta,pi,deg2rad
  double complex R(3,3)
  double precision c,s,ux,uy,uz,norme_axe
  integer i,j

  pi=3.1415926535897932d0
  deg2rad = pi / 180.d0
  ! on normalise l'axe
  !write(100,*) 'axe rotation  :',axe
  norme_axe=dsqrt(axe(1)**2 + axe(2)**2 + axe(3)**2)

  ! composantes de l'axe
  ux=axe(1)/norme_axe
  uy=axe(2)/norme_axe
  uz=axe(3)/norme_axe

  ! on calcule le cos et sin
  c=dcos(deg2rad * theta);s=dsin(deg2rad * theta)

  ! matrice de rotation complexe
  R(1,1)=dcmplx(ux**2 + (1.d0-ux**2)*c)
  R(1,2)=dcmplx(ux*uy*(1.d0-c)-uz*s)
  R(1,3)=dcmplx(ux*uz*(1.d0-c)+uy*s)

  R(2,1)=dcmplx(ux*uy*(1.d0-c)+uz*s)
  R(2,2)=dcmplx(uy**2+(1.d0-uy**2)*c)
  R(2,3)=dcmplx(uy*uz*(1.d0-c)-ux*s)

  R(3,1)=dcmplx(ux*uz*(1.d0-c)-uy*s)
  R(3,2)=dcmplx(uy*uz*(1.d0-c)+ux*s)
  R(3,3)=dcmplx(uz**2+(1.d0-uz**2)*c)


end subroutine rotation_matrix

!-------------------------------------------------------------------------------
! R=R2*R1*R0
subroutine compose3matrix(R,R0,R1,R2)
  implicit none
  double complex R(3,3),R0(3,3),R1(3,3),R2(3,3)
  integer i,j,k

  R(:,:)=dcmplx(0.d0)
 ! multiplication R=R1*R0
  do j=1,3
   do i=1,3
    do k=1,3
       R(i,j)=R(i,j) + R1(i,k)*R0(k,j)
    enddo
   enddo
  enddo

  R1(:,:)=R(:,:)
  R(:,:)=dcmplx(0.d0)
! multiplication R=R2*R1
  do j=1,3
   do i=1,3
    do k=1,3
        R(i,j)=R(i,j) + R2(i,k)*R1(k,j)
    enddo
   enddo
  enddo
end subroutine compose3matrix

!------------------------------------------------------------------------------
! rotation pour passer d'un repere local a un autre
subroutine local2localMatrix(lat0,lon0,lat1,lon1,R)
  implicit none
  double precision lat0,lon0,lat1,lon1
  double precision distance_epicentrale,azi,bazi
  double complex R(3,3),axe_rotation(3)
  ! calcul de la distance epicentrale = angle de rotation
  call epitra1(lat0,lon0,lat1,lon1,distance_epicentrale,azi,bazi)
  ! calcul de l'axe de rotation = perendiculaire au plan (O,P0,P1)
  call calcule_axe_rotation(lat0,lon0,lat1,lon1,axe_rotation)
  ! on calcule la matrice de rotation
  call rotation_matrix(R,axe_rotation,distance_epicentrale)
end subroutine local2localMatrix

!-------------------------------------------------------------------------------
!calcul de l'axe de rotation
subroutine calcule_axe_rotation(lat0,lon0,lat1,lon1,axe)
  implicit none
  double precision lat0,lon0,lat1,lon1,axe(3)
  double precision X0(3),X1(3)
  ! on passe dans le repere global cartesien
  call geograph2cartglob(X0,lat0,lon0,1.d0)
  call geograph2cartglob(X1,lat1,lon1,1.d0)
  ! on fait le produit vectoriel X0^X1
  call pdt_vectoriel(axe,X0,X1)
end subroutine calcule_axe_rotation

!-------------------------------------------------------------------------------
! passage geographique -> global
subroutine geograph2cartglob(X,lat,lon,r)
  implicit none
  double precision deg2rad,lat,lon,r,X(3)
  integer i
  deg2rad=3.1415926535897932d0/180.d0
  X(1)=r*dcos(deg2rad*lon)*cos(deg2rad*lat);
  X(2)=r*dsin(deg2rad*lon)*cos(deg2rad*lat);
  X(3)=r*dsin(deg2rad*lat);
end subroutine geograph2cartglob

!-------------------------------------------------------------------------------
! passage global -> geographique
subroutine cartglob2geograph(X,lat,lon,r)
  implicit none
  double precision r,lat,lon,X(3),rad2deg
  rad2deg=180.d0/3.1415926535897932d0
  r=dsqrt(X(1)**2+X(2)**2+X(3)**2);
  lon=datan2(X(2),X(1))*rad2deg;
  lat=dasin(X(3)/r)*rad2deg;
end subroutine cartglob2geograph

!-------------------------------------------------------------------------------
! produit vectoriel
subroutine pdt_vectoriel(Z,X,Y)
  implicit none
  double precision X(3),Y(3),Z(3)
  z(1)=x(2)*y(3)-x(3)*y(2);
  z(2)=x(3)*y(1)-x(1)*y(3);
  z(3)=x(1)*y(2)-x(2)*y(1);
end subroutine pdt_vectoriel

!-------------------------------------------------------------------------------
! produit matrice vecteur Y=R*X
subroutine matmulvect(Y,R,X)
  implicit none
  double complex Y(3),R(3,3),X(3)
  integer i,k
  Y(:)=dcmplx(0.d0)
  do i=1,3
    do k=1,3
      Y(i)=Y(i)+R(i,k)*X(k)
    enddo
  enddo
end subroutine matmulvect
!------------------------------------------------------------------------------
! affichage d'une matrice complexe
subroutine Display_matrix_complex(iunit,M)
  implicit none
  integer i,iunit
  double complex M(3,3)
  do i=1,3
  write(iunit,'("|",2f10.5,5x,2f10.5,5x,2f10.5," |")') real(M(i,1)),aimag(M(i,1)),real(M(i,2)),aimag(M(i,2)),real(M(i,3)),aimag(M(i,3))
  enddo
end subroutine Display_matrix_complex
!------------------------------------------------------------------------------
! affichage d'une matrice complexe partie reele
subroutine Display_matrix_realpart(iunit,M)
  implicit none
  integer i,iunit
  double complex M(3,3)
  do i=1,3
  write(iunit,'("|",f10.5,5x,f10.5,5x,f10.5," |")') real(M(i,1)),real(M(i,2)),real(M(i,3))
  enddo
end subroutine Display_matrix_realpart

