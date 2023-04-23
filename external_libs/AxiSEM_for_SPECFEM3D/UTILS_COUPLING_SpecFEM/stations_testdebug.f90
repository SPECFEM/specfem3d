program stations_testdebug

  implicit none

  double precision  ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
  double precision  lon_center_chunk, lat_center_chunk, chunk_azi
  double precision  chunk_depth,deg2rad
  double precision  rotation_matrix(3,3),rotation_inv_matrix(3,3)
  character(len=100) line
  character(len=5) name_sta
  character(len=1) code

  double precision x, y, z, z_bot
  double precision lat, long, radius
  integer k

  code = 'S'

  deg2rad = 3.141592653589793d0/180.d0

!

  open(10, file='ParFileMeshChunk')

  read(10,'(a)') line
  read(10,*) ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
  read(10,'(a)') line
  read(10,*) lon_center_chunk, lat_center_chunk, chunk_azi
  read(10,'(a)') line
  read(10,*) chunk_depth

  close(10)

!

  ANGULAR_WIDTH_XI_RAD  = deg2rad * ANGULAR_WIDTH_XI_RAD
  ANGULAR_WIDTH_ETA_RAD = deg2rad * ANGULAR_WIDTH_ETA_RAD
  chunk_depth           = chunk_depth * 1000.d0

  call compute_rotation_matrix(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)
  call compute_inv_rotation_matrix(rotation_inv_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)

  z_bot = 6120505.36002571 !! to be found in file MESH/model_1D.in

  open(11,file='input_testbox_to_convert_CART.txt')
  open(12,file='input_testbox_to_convert_SPHER.txt')
  open(13,file='input_testbox_to_convert_RECART.txt')


  k = 0

  do

     k = k + 1

!     write(name_sta,'(a,i4.4)') code, k

     read(11,*,end=99) x, y, z

     call cart2geogr(x, y, z, rotation_matrix, long, lat, radius)

     write(*,*) radius, lat, long
     write(12,'(4f20.10)') radius, lat, long


     call geogr2cart(x, y, z, rotation_inv_matrix, long, lat, radius)

     write(*,*) x, y, z
     write(13,*) x, y, z

  enddo


99 close(11)

  close(12)
  close(13)

  !lat=1.d0
  !long=61.d0
  !radius=6371.d0
  !call geogr2cart(x,y,z,rotation_inv_matrix,long,lat,radius)
  !write(*,*) x,y,z

end program stations_testdebug

!=======================================================================================================

  subroutine cart2geogr(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius)

    implicit none

    integer NDIM,NGLLX,NGLLY,NGLLZ
    double precision xstore,ystore,zstore
    double precision longitud,latitud,radius
    double precision rotation_matrix(3,3)
    double precision vector_ori(3),vector_rotated(3)
    double precision rayon,x,y,z,deg2rad,long,lati
    integer i,j,k

    deg2rad=3.141592653589793d0/180.d0

    NDIM=3

    vector_ori(1)=xstore
    vector_ori(2)=ystore
    vector_ori(3)=zstore

    do i = 1,NDIM  !
       vector_rotated(i) = 0.d0
       do j = 1,NDIM

          vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)

       enddo
    enddo

    x        = vector_rotated(1)
    y        = vector_rotated(2)
    z        = vector_rotated(3)
    rayon    = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

    long     = atan2(y,x)
    lati     = asin(z/rayon)

    longitud = long/deg2rad
    latitud  = lati/deg2rad
    radius   = rayon/1000.d0

  end subroutine cart2geogr


!=======================================================================================================



!=======================================================================================================
!!!
!!
  subroutine geogr2cart(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius)
!
    implicit none

    integer NDIM
    double precision xstore,ystore,zstore
    double precision longitud,latitud,radius
    double precision rotation_matrix(3,3)
    double precision vector_ori(3),vector_rotated(3)
    double precision rayon,x,y,z,deg2rad,long,lati
    integer i,j,k

    deg2rad=3.141592653589793d0/180.d0
    NDIM=3

    vector_ori(1)=1000.d0*radius*cos(deg2rad*longitud)*sin(deg2rad*(90.d0-latitud))
    vector_ori(2)=1000.d0*radius*sin(deg2rad*longitud)*sin(deg2rad*(90.d0-latitud))
    vector_ori(3)=1000.d0*radius*cos(deg2rad*(90.d0-latitud))

    write(*,*) 'v3', vector_ori(3)

    do i = 1,NDIM
       vector_rotated(i) = 0.d0
       do j = 1,NDIM

          vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)

       enddo
    enddo

    xstore=vector_rotated(1)
    ystore=vector_rotated(2)
    zstore=vector_rotated(3)

  end subroutine geogr2cart



!=======================================================================================================
!
!=======================================================================================================

  subroutine  compute_rotation_matrix(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)

  implicit none

  double precision rotation_matrix(3,3),lon_center_chunk,lat_center_chunk, chunk_azi
  double precision R0(3,3),R1(3,3),R2(3,3),axe_rotation(3),R00(3,3)

  ! je met le chunk en 0,0
  axe_rotation(1)=0.d0; axe_rotation(2)=1.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R00,axe_rotation,90.d0)  ! je ramene le chunk en (0,0)
  ! rotation de l'azimuth du chunk
  axe_rotation(1)=1.d0; axe_rotation(2)=0.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R0,axe_rotation,90.-chunk_azi)
  ! on met le chunk a la bonne latitude
  axe_rotation(1)=0.d0; axe_rotation(2)=-1.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R1,axe_rotation,lat_center_chunk)
  ! on met le chunk a la bonne longitude
  axe_rotation(1)=0.d0; axe_rotation(2)=0.d0; axe_rotation(3)=1.d0
  call rotation_matrix_axe(R2,axe_rotation, lon_center_chunk)
  ! rotation resultante
  call compose4matrix(rotation_matrix,R00,R0,R1,R2)

  end subroutine compute_rotation_matrix

!=======================================================================================================

  subroutine  compute_inv_rotation_matrix(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)

  implicit none

  double precision rotation_matrix(3,3),lon_center_chunk,lat_center_chunk, chunk_azi
  double precision R0(3,3),R1(3,3),R2(3,3),axe_rotation(3),R00(3,3)

  ! on met le chunk a la bonne longitude
  axe_rotation(1)=0.d0; axe_rotation(2)=0.d0; axe_rotation(3)=-1.d0
  call rotation_matrix_axe(R00,axe_rotation, lon_center_chunk)
  ! on met le chunk a la bonne latitude
  axe_rotation(1)=0.d0; axe_rotation(2)=1.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R0,axe_rotation,lat_center_chunk)
  ! rotation de l'azimuth du chunk
  axe_rotation(1)=-1.d0; axe_rotation(2)=0.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R1,axe_rotation,90.-chunk_azi)
! je met le chunk en 0,0
  axe_rotation(1)=0.d0; axe_rotation(2)=-1.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R2,axe_rotation,90.d0)  ! je ramene le chunk en (0,0)
  ! rotation resultante
  call compose4matrix(rotation_matrix,R00,R0,R1,R2)

end subroutine compute_inv_rotation_matrix


!=======================================================================================================
!
!   ROUTINES POUR FAIRE DES ROTATIONS 3D ET DIVERS CHANGEMENTS DE REPERES
!
! Vadim Monteiller Mars 2013
!
!-------------------------------------------------------------------------------
! matrice de rotation 3D d'axe "axe" et d'angle theta (en degres)
! cette matrice est en complexe
!
!=======================================================================================================
!
  subroutine rotation_matrix_axe(R,axe,theta)

  implicit none

  double precision axe(3),theta,pi,deg2rad
  double precision R(3,3)
  double precision c,s,ux,uy,uz,norme_axe

  pi=3.1415926535897932d0
  deg2rad = pi / 180.d0
  ! on normalise l'axe
  norme_axe=dsqrt(axe(1)**2 + axe(2)**2 + axe(3)**2)

  ! composantes de l'axe
  ux=axe(1)/norme_axe
  uy=axe(2)/norme_axe
  uz=axe(3)/norme_axe

  ! on calcule le cos et sin
  c=dcos(deg2rad * theta);s=dsin(deg2rad * theta)

  ! matrice de rotation complexe
  R(1,1)=(ux**2 + (1.d0-ux**2)*c)
  R(1,2)=(ux*uy*(1.d0-c)-uz*s)
  R(1,3)=(ux*uy*(1.d0-c)+uy*s)

  R(2,1)=(ux*uy*(1.d0-c)+uz*s)
  R(2,2)=(uy**2+(1.d0-uy**2)*c)
  R(2,3)=(uy*uz*(1.d0-c)-ux*s)

  R(3,1)=(ux*uz*(1.d0-c)-uy*s)
  R(3,2)=(uy*uz*(1.d0-c)+ux*s)
  R(3,3)=(uz**2+(1.d0-uz**2)*c)

  !write(49,*) ' MATRICE ROTATION '
  !write(49,*) R(1,:)
  !write(49,*) R(2,:)
  !write(49,*) R(3,:)
  !write(49,*)

  end subroutine rotation_matrix_axe

!=======================================================================================================

!=======================================================================================================
!
! R=R2*R1*R0
!
!=======================================================================================================

  subroutine compose4matrix(R,R00,R0,R1,R2)

  implicit none

  double precision R(3,3),R0(3,3),R1(3,3),R2(3,3),R00(3,3),Rtmp(3,3)
  integer i,j,k


  R(:,:)=0.d0
  ! multiplication R=R0*R00
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R0(i,k)*R00(k,j)
        enddo
     enddo
  enddo

  ! multiplication R=R1*R
  Rtmp=R
  R(:,:)=0.d0
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R1(i,k)*Rtmp(k,j)
        enddo
     enddo
  enddo

  ! multiplication R=R2*R
  Rtmp=R
  R(:,:)=0.d0
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R2(i,k)*Rtmp(k,j)
        enddo
     enddo
  enddo

  !write(49,*) ' MATRICE ROTATION COMPLETE '
  !write(49,*) R(1,:)
  !write(49,*) R(2,:)
  !write(49,*) R(3,:)
  !write(49,*)

  end subroutine compose4matrix
