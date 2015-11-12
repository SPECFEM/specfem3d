program rotation_specfem_solution

  implicit none
  double precision  ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
  double precision  lon_center_chunk, lat_center_chunk, chunk_azi
  double precision  chunk_depth,deg2rad
  double precision  rot_matrix(3,3),rot_inv_matrix(3,3)
  double precision            :: rot(3,3), irot(3,3)
  double precision            :: rot_chunk(3,3), irot_chunk(3,3)
  character(len=100) line
  character(len=5) name_sta
  character(len=1) code
  character(len=120) sis_cx,sis_cy,sis_cz,output_file_E, output_file_N, output_file_Z
  double precision x,y,z,z_bot
  double precision t,lat,long,radius
  integer k

  code='S'

  deg2rad = 3.141592653589793d0/180.d0
  open(10, file='ParFileMeshChunk')

  read(10,'(a)') line
  read(10,*) ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
  read(10,'(a)') line
  read(10,*) lon_center_chunk, lat_center_chunk, chunk_azi
  read(10,'(a)') line
  read(10,*) chunk_depth

  close(10)

  ANGULAR_WIDTH_XI_RAD  = deg2rad * ANGULAR_WIDTH_XI_RAD
  ANGULAR_WIDTH_ETA_RAD = deg2rad * ANGULAR_WIDTH_ETA_RAD
  chunk_depth           = chunk_depth * 1000.d0

!!$  call compute_rotation_matrix(rotation_matrix, lon_center_chunk, lat_center_chunk, chunk_azi)
!!$  call compute_inv_rotation_matrix(rotation_inv_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)

  !!! WARNING :: chunk_azi is not taken into account  !!!!!!!!!!!
  call rotation_matrix(lon_center_chunk, lat_center_chunk, rot_chunk, irot_chunk)
!!$  write(*,*) rotation_matrix(1,:)
!!$  write(*,*) rotation_matrix(2,:)
!!$  write(*,*) rotation_matrix(3,:)
!!$
!!$  write(*,*)
!!$
!!$
!!$  write(*,*) rotation_inv_matrix(1,:)
!!$  write(*,*) rotation_inv_matrix(2,:)
!!$  write(*,*) rotation_inv_matrix(3,:)
!!$
!!$  write(*,*)
!!$
!!$  lat=0.d0
!!$  long=60.d0
!!$  radius=6371.d0
  read(*,*) z_bot  !! to be found in file MESH/model_1D.in
  open(10,file='stations_to_convert.txt')
  !open(21,file='station_converted.txt')
  read(*,'(a)') sis_cx
  read(*,'(a)') sis_cy
  read(*,'(a)') sis_cz
  read(*,'(a)') output_file_E
  read(*,'(a)') output_file_N
  read(*,'(a)') output_file_Z
  open(21,file=trim(sis_cx))
  open(22,file=trim(sis_cy))
  open(23,file=trim(sis_cz))
  open(24,file=trim(output_file_E))
  open(25,file=trim(output_file_N))
  open(26,file=trim(output_file_Z))
  k=0
  !do
     k=k+1
     write(name_sta,'(a,i4.4)') code,k
     read(10,*,end=99) radius,lat,long
     call rotation_matrix(long, lat, rot, irot)

     call product_matrix(irot,rot_chunk, rot_matrix)
     do
        read(21,*,end=98) t,x
        read(22,*,end=98) t,y
        read(23,*,end=98) t,z
        call multip_mv(rot_matrix, x,y,z, long, lat, radius)
        !call cart2geogr(x,y,z,rotation_matrix,long,lat,radius)
        write(24,'(2f30.15)') t + 703.80324d0, long
        write(25,'(2f30.15)') t + 703.80324d0, lat
        write(26,'(2f30.15)') t + 703.80324d0, radius
     enddo
     98 continue
     close(21)
     close(22)
     close(23)
     close(24)
     close(25)
     close(26)
     !write(21,'(a5,1x,a2,1x,4f20.5)') name_sta,'SY',y,x,-(z-z_bot),-(z-z_bot) !! because specfem is reading y and x actually
  !enddo
99 close(10)
  !close(21)
  !lat=1.d0
  !long=61.d0
  !radius=6371.d0

  !call geogr2cart(x,y,z,rotation_inv_matrix,long,lat,radius)
  !write(*,*) x,y,z


end program rotation_specfem_solution



!=======================================================================================================
!!!
!!
  subroutine cart2geogr(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius)
!
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

    x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
    rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

    long=atan2(y,x)
    lati=asin(z/rayon)

    longitud = long/deg2rad
    latitud  = lati/deg2rad
    radius   = rayon/1000.d0




  end subroutine cart2geogr


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

  write(49,*) ' MATRICE ROTATION '
  write(49,*) R(1,:)
  write(49,*) R(2,:)
  write(49,*) R(3,:)
  write(49,*)

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

  write(49,*) ' MATRICE ROTATION COMPLETE '
  write(49,*) R(1,:)
  write(49,*) R(2,:)
  write(49,*) R(3,:)
  write(49,*)

  end subroutine compose4matrix


  subroutine rotation_matrix(long, lat, rot, irot)

    double precision            :: long, lat
    double precision            :: rot(3,3), irot(3,3)
    double precision            :: ct, st, cp, sp
    double precision, parameter :: deg2rad=0.017453292519943d0

    ct = cos(deg2rad * long)
    st = sin(deg2rad * long)

    cp = cos(deg2rad * (90.d0 - lat))
    sp = sin(deg2rad * (90.d0 - lat))


    !! rotation matrix
    rot(1,1) = -st ; rot(1,2) = -cp*ct; rot(1,3) = ct*sp
    rot(2,1) = ct  ; rot(2,2) = -cp*st; rot(2,3) = sp*st
    rot(3,1) = 0.d0; rot(3,2) = sp    ; rot(3,3) = cp

    !! inverse=transpose
    irot(:,:)=rot(:,:)
    irot(1,2)=rot(2,1);irot(1,3)=rot(3,1)
    irot(2,1)=rot(1,2);irot(2,3)=rot(3,2)
    irot(3,1)=rot(1,3);irot(3,2)=rot(2,3)

  end subroutine rotation_matrix

  subroutine  product_matrix(irot,rot_chunk, rotation_matrix)

    integer     :: i,j,k
    double precision            :: irot(3,3), rotation_matrix(3,3),rot_chunk(3,3)

    rotation_matrix(:,:)=0.d0

    do i=1,3
       do j=1,3
          do k=1,3
             rotation_matrix(i,j)= rotation_matrix(i,j)+irot(i,k)*rot_chunk(k,j)
          enddo
       enddo
    enddo


  end subroutine product_matrix

  subroutine multip_mv(rotation_matrix, x,y,z, long, lat, radius)

    double precision     :: x, y, z
    double precision     :: long, lat, radius
    double precision     :: rotation_matrix(3,3)
    double precision     :: vect(3), vect1(3)

    vect1(1)=x
    vect1(2)=y
    vect1(3)=z

    vect(:)=0.d0
    do i=1,3
       do j=1,3
          vect(i)=vect(i)+rotation_matrix(i,j)*vect1(j)
       enddo
    enddo

    long=vect(1)
    lat=vect(2)
    radius=vect(3)

  end subroutine multip_mv
