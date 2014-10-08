!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine chunk_of_earth_Mesh()

  implicit none

!==============================================================================================!
!                                                                                              !
!  Singular option of meshfem3D : MESH OF A GLOBE EARTH CHUNK FOR THE INTERFACE DSM-SPECFEM3D  !
!  Case of 8 nodes per element (HEX8)                                                          !
!                                                                                              !
!  Vadim Monteiller, February 2013                                                             !
!  Integrated in meshfem3d by CD, September 2014                                               !
!                                                                                              !
!  WARNING : A local convention is used for the mapping of                                     !
!            the cubic sphere (to complete)                                                    !
!                                                                                              !
!==============================================================================================!

!
!--- Parameters
!

  integer, parameter :: NGLLX  = 5, NGLLY = 5, NGLLZ = 5, NGNOD = 8, NDIM = 3
  integer, parameter :: myrank = 0
  integer, parameter :: nlayer = 12 !! (number of layer in the model iasp91, or ak135, or prem (one more layer than the model)

  double precision, parameter :: GAUSSALPHA = 0.d0, GAUSSBETA = 0.d0

  logical, parameter ::  RUN_BENCHMARK = .false.

!
!--- Other
!

  integer  nel_lat, nel_lon, nel_depth, NX, NY, NZ, Ndepth, nglob, kglob, ilocnum, ieoff, npointot
  integer ilat, ilon, ispec, iz, i, j, k, nspec, ia, izshift, index_mat
  integer ispec2Dxmin, ispec2Dxmax, ispec2Dymin, ispec2Dymax, ispec2Dzmin, ispec2Dzmax
  integer ilayer_current, ilayer
  integer nlat_dsm, nlon_dsm

  integer iaddx(NGNOD), iaddy(NGNOD), iaddz(NGNOD)

  integer, allocatable :: inum_loc(:,:,:,:), iglob(:), loc(:), current_layer(:)

  double precision ratio_eta, ratio_xi
  double precision ANGULAR_WIDTH_ETA_RAD, ANGULAR_WIDTH_XI_RAD, Z_DEPTH_BLOCK, UTM_X_MIN, UTM_X_MAX
  double precision lat_center_chunk, lon_center_chunk, chunk_depth, chunk_azi
  double precision R_EARTH, TINYVAL, PI, ZERO, deg2rad
  double precision x, y, z, px, py, pz, z_bottom

  double precision rotation_matrix(3,3)
  double precision zlayer(nlayer), vpv(nlayer,4), vsv(nlayer,4), density(nlayer,4)
  double precision xelm(NGNOD), yelm(NGNOD), zelm(NGNOD)
  double precision xstore(NGLLX,NGLLY,NGLLZ), ystore(NGLLX,NGLLY,NGLLZ), zstore(NGLLX,NGLLY,NGLLZ)

  !! 3D shape functions and their derivatives
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ), dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  !! GLL points and weights of integration
  double precision xigll(NGLLX), yigll(NGLLY), zigll(NGLLZ), wxgll(NGLLX), wygll(NGLLY), wzgll(NGLLZ)

  double precision, allocatable :: xp(:), yp(:), zp(:), xgrid(:,:,:,:), ygrid(:,:,:,:), zgrid(:,:,:,:)
  double precision, allocatable :: lon_zmin(:,:), lat_zmin(:,:)
  double precision, dimension(:,:), allocatable :: ProfForGemini

  logical test

  logical, allocatable :: ifseg(:)
  logical, dimension(:,:), allocatable :: iboun ! boundary locator

  character(len=100) line
  character(len=250) model1D_file

  character(len=10), parameter :: MESH = "./MESH/"

!! Unused
!
! integer iii, jjj, kkk,
! double precision long, lati, x_bot, y_bot, z_bot,rayon,ratio_depth,
! double precision theta, colat, vector_ori(3), vector_rotated(3)

!
!--- WARNING ==> CONVENTION : (lon,lat) -> (xi,eta)
!---                          (k = 6 with -z for the mapping of the cubic sphere, cf Chervot 2012)
!---                          We define the mesh of a chunk of the earth in the cubic sphere
!

  PI      = 3.141592653589793d0
  deg2rad = 3.141592653589793d0/180.d0
  R_EARTH = 6371000.d0
  TINYVAL = 1.d-9
  ZERO    = 0.d0

  open(49, file=trim(MESH)//'output_mesher_chunk.txt')

  if (RUN_BENCHMARK) then
!
!--- Parameters fixed inside the subroutine for the moment
!
     ANGULAR_WIDTH_ETA_RAD = 10.d0 * deg2rad ! latitude 2.5
     ANGULAR_WIDTH_XI_RAD  = 20.d0 * deg2rad ! longitude 2.0

!    Chunk center
     lat_center_chunk      = 0.d0  ! 42.35d0 ! 42.5d0 !* deg2rad
     lon_center_chunk      = 60.d0 ! 1.3d0   ! 1.2d0  !* deg2rad

!    Azimuth
     chunk_azi             = 0.d0 !90.d0 !80.d0 !10.d0  !* deg2rad

!    Depth
     chunk_depth           = 1000.d0 * 1000.d0 ! 250.d0 * 1000.d0

!    Number of elements
     nel_lat               = 20 ! 120  ! 15
     nel_lon               = 40 ! 96   ! 15
     nel_depth             = 20 ! 100  ! 10

  else

     open(10, file=trim(MESH)//'ParFileMeshChunk')

     read(10,'(a)') line
     read(10,*) ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
     read(10,'(a)') line
     read(10,*) lon_center_chunk, lat_center_chunk, chunk_azi
     read(10,'(a)') line
     read(10,*) chunk_depth
     read(10,'(a)') line
     read(10,*) nel_lon,nel_lat, nel_depth
     read(10,'(a)') line
     read(10,'(a)') model1D_file

     model1D_file = 'MESH/'//trim(model1D_file)

     close(10)

     ANGULAR_WIDTH_XI_RAD  = deg2rad * ANGULAR_WIDTH_XI_RAD
     ANGULAR_WIDTH_ETA_RAD = deg2rad * ANGULAR_WIDTH_ETA_RAD
     chunk_depth           = chunk_depth * 1000.d0

  endif

  NX = nel_lon
  NY = nel_lat
  NZ = nel_depth

!
!===========================================================================
!
!--- TO DO : the reference chunk must be always symmetric (EW) and (NS)
!

  nlon_dsm = (ngllx - 1) * NX + 1
  nlat_dsm = (nglly - 1) * NY + 1
  nglob    = (nel_lat + 1) * (nel_lon + 1) * (nel_depth + 1)
  nspec    = nel_lat * nel_lon * nel_depth
  npointot = 8 * nspec

  allocate(xp(npointot), yp(npointot), zp(npointot))
  allocate(iglob(npointot), loc(npointot))
  allocate(ifseg(npointot))
  allocate(ProfForGemini(0:NZ-1,3))
  allocate(current_layer(0:NZ-1))
  allocate(inum_loc(2,2,2,nspec))
  allocate(xgrid(2,2,2,nspec), ygrid(2,2,2,nspec), zgrid(2,2,2,nspec))
  allocate(lon_zmin(nlon_dsm,nlat_dsm), lat_zmin(nlon_dsm,nlat_dsm))
  allocate(iboun(6,nspec)) ! boundary locator

  iboun(:,:) = .false.

  iaddx(1) = 0
  iaddy(1) = 0
  iaddz(1) = 0

  iaddx(2) = 1
  iaddy(2) = 0
  iaddz(2) = 0

  iaddx(3) = 1
  iaddy(3) = 1
  iaddz(3) = 0

  iaddx(4) = 0
  iaddy(4) = 1
  iaddz(4) = 0

  iaddx(5) = 0
  iaddy(5) = 0
  iaddz(5) = 1

  iaddx(6) = 1
  iaddy(6) = 0
  iaddz(6) = 1

  iaddx(7) = 1
  iaddy(7) = 1
  iaddz(7) = 1

  iaddx(8) = 0
  iaddy(8) = 1
  iaddz(8) = 1

!
!===========================================================================
!
!--- set up coordinates of the Gauss-Lobatto-Legendre points
!
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

!
!--- if number of points is odd, the middle abscissa is exactly zero
!
  if(mod(NGLLX,2) /= 0) xigll((NGLLX - 1)/2 + 1) = ZERO
  if(mod(NGLLY,2) /= 0) yigll((NGLLY - 1)/2 + 1) = ZERO
  if(mod(NGLLZ,2) /= 0) zigll((NGLLZ - 1)/2 + 1) = ZERO

!
!--- get the 3-D shape functions
!
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll,NGNOD)

!
!--- rotation matrix to switch to the geographical coordinates
!--- call euler_angles(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)
!--- new rotation matrix
!

  call compute_rotation_matrix(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)

!
!--- call ReadIasp91(vpv,vsv,density,zlayer,nlayer)
!

  call Read_dsm_model(model1D_file,vpv,vsv,density,zlayer,nlayer)

!
!--- calculation of the vertical discretization of layers
!

  Z_DEPTH_BLOCK = chunk_depth/1000.d0 !!!! switch to km

  call CalGridProf(ProfForGemini,current_layer,zlayer,nlayer,NZ,Z_DEPTH_BLOCK)

!
!===========================================================================
!
!--- GRID OF THE MESH
!
  izshift     = 0
  ispec       = 0
  kglob       = 0
  Ndepth      = 0
  ispec2Dxmin = 0
  ispec2Dxmax = 0
  ispec2Dymin = 0
  ispec2Dymax = 0
  ispec2Dzmin = 0
  ispec2Dzmax = 0

!
!--- Interface file  DSM-SPECFEM3D
!

  open(27, file = trim(MESH)//'.recdepth')  ! receptors on the vertical
  open(28, file = trim(MESH)//'stxmin')
  write(28,*) nlat_dsm          ! face xmin

  open(29, file = trim(MESH)//'stxmax')
  write(29,*) nlat_dsm ! face xmax

  open(30, file = trim(MESH)//'stymin')
  write(30,*) nlon_dsm ! face ymin

  open(31, file = trim(MESH)//'stymax')
  write(31,*) nlon_dsm ! face ymax

  open(38, file = trim(MESH)//'IgXmin')
  open(39, file = trim(MESH)//'IgXmax')
  open(40, file = trim(MESH)//'IgYmin')
  open(41, file = trim(MESH)//'IgYmax')
  open(42, file = trim(MESH)//'IgZmin')

! MESH for SPECFEM3D
  open(86, file = trim(MESH)//'nummaterial_velocity_file')
  open(87, file = trim(MESH)//'materials_file')

! open(88, file = 'model_1D.in')
  open(88, file = trim(MESH)//'OrigRepSpecfm')
  write(88,*) lon_center_chunk, lat_center_chunk
  write(88,*) chunk_azi, ANGULAR_WIDTH_XI_RAD/deg2rad, ANGULAR_WIDTH_ETA_RAD/deg2rad
  close(88)

  open(89, file = trim(MESH)//'flags_boundary.txt')
  open(90, file = trim(MESH)//'Nb_ielm_faces.txt')

! open(32,file='gll_zmin')
! open(125,file='ggl_elemts')

!
!-- Loop on the grid of the spectral elements
!

  ilayer    = 0
  index_mat = 0

  do iz = 0, nel_depth - 1

     ilayer_current=current_layer(iz) - 1 ! Caution between pickets and intervals !!

     if (iz /= 0) then
        if (current_layer(iz-1) /= current_layer(iz)) then
           izshift   = izshift + 1 ! point is repeated on the interface for DSM
           index_mat = index_mat - 1
           write(86,'(a1,2x,i10,2x,a10,2x,a7,2x,a20,2x,a1)') &
                '2', index_mat, 'tomography', 'elastic', 'tomography_model.xyz', '1'
        endif

     else
!       We write the first material
        index_mat = index_mat - 1
        write(86,'(a1,2x,i10,2x,a10,2x,a7,2x,a20,2x,a1)') &
             '2', index_mat, 'tomography', 'elastic', 'tomography_model.xyz', '1'
     endif

     do ilat=0,nel_lat-1
        do ilon=0,nel_lon-1

           ispec = ispec + 1
           ! material file
           write(87 ,*) ispec,index_mat

           ! get boundary
           ! on boundary 1: x=xmin
           if(ilon == 0 ) then
              iboun(1,ispec)=.true.
              ispec2Dxmin=ispec2Dxmin+1
              write(89,*) ispec,ispec2Dxmin,1
           endif
           ! on boundary 2: xmax
           if(ilon == nel_lon-1) then
              iboun(2,ispec)=.true.
              ispec2Dxmax=ispec2Dxmax+1
              !write(*,*) '------ TOZ',ispec,ilon
              write(89,*) ispec,ispec2Dxmax,2
           endif
           ! on boundary 3: ymin
           if(ilat == 0) then
              iboun(3,ispec)=.true.
              ispec2Dymin=ispec2Dymin+1
              write(89,*) ispec,ispec2Dymin,3
           endif
           ! on boundary 4: ymax
           if(ilat == nel_lat-1 ) then
              iboun(4,ispec) =.true.
              ispec2Dymax=ispec2Dymax+1
              write(89,*) ispec,ispec2Dymax,4
           endif
           ! on boundary 5: bottom
           if(iz == 0) then
              iboun(5,ispec)=.true.
              ispec2Dzmin=ispec2Dzmin+1
              write(89,*) ispec,ispec2Dzmin,5
           endif
           ! on boundary 6: top
           if(iz == nel_depth-1) then
              ispec2Dzmax= ispec2Dzmax+1
              iboun(6,ispec)=.true.
           endif

          ! 8 vertices of the element ispec
           do ia=1,NGNOD

              i=iaddx(ia)
              j=iaddy(ia)
              k=iaddz(ia)

              z = 1000d0*ProfForGemini(iz,1+k)

              ! longitude
              ratio_xi = (dble(ilon+i)) / dble(NX)
              x = 2.d0*ratio_xi-1.d0
              x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

              ! latitude
              ratio_eta = (dble(ilat+j)) / dble(NY)
              y = 2.d0*ratio_eta-1.d0
              y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

              !if (ilat==0.and.iz==0) write(49,*) ia,i,ratio_xi
              !mapping cubic sphere (k=5) (Chevrot et al 2012)
              ! (opposite sign)
              !pz = z/dsqrt(1.d0 + y*y + x*x)
              !px = x*pz
              !py = y*pz
              ! mapping which makes a chunk at the North Pole
              ! mapping cubic sphere (k=6, Chevrot at al 2012, avec -z)

              pz=  z/dsqrt(1.d0 + y*y + x*x) !(=r/s)
              px= pz * x !(tan(xi) * r/s)
              py= pz * y !(tan(eta) * r/s)

              ! old version
              xgrid(i+1,j+1,k+1,ispec) = px !px
              ygrid(i+1,j+1,k+1,ispec) = py !py
              zgrid(i+1,j+1,k+1,ispec) = pz

              !xgrid(i+1,j+1,k+1,ispec) = py ! long
              !ygrid(i+1,j+1,k+1,ispec) = pz ! lat
              !zgrid(i+1,j+1,k+1,ispec) = px ! prof

              xelm(ia)=xgrid(i+1,j+1,k+1,ispec)
              yelm(ia)=ygrid(i+1,j+1,k+1,ispec)
              zelm(ia)=zgrid(i+1,j+1,k+1,ispec)


           enddo

           ! INTERFACE FOR DSM ------

           ! Vertical receptors
           if (ilat==0 .and. ilon==0) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              call write_gllz_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,current_layer,nel_depth,ilayer,iz,Ndepth)
           endif

           ! Horizontal receptors

           ! stxmin
           if (ilon==0.and.iz==nel_depth-1)  then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              if (ilat==nel_lat-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
              call  write_stxmin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

           endif
            if (ilon==0) call write_Igm_file(38,ispec2Dxmin,NGLLY,NGLLZ,ilat,iz,izshift,ilayer_current)

           ! stxmax
           if (ilon==nel_lon - 1 .and. iz==nel_depth-1)  then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
               if (ilat==nel_lat-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
              call  write_stxmax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)
           endif
           if (ilon==nel_lon-1)  call write_Igm_file(39,ispec2Dxmax,NGLLY,NGLLZ,ilat,iz,izshift,ilayer_current)

           ! stymin
           if (ilat==0.and. iz==nel_depth-1)  then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              if (ilon==nel_lon-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
               call write_stymin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)
           endif
           if (ilat==0) call write_Igm_file(40,ispec2Dymin,NGLLX,NGLLZ,ilon,iz,izshift,ilayer_current)
           ! stymax
           if (ilat==nel_lat-1.and. iz==nel_depth-1)  then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              if (ilon==nel_lon-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
              call write_stymax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)
           endif
           if (ilat==nel_lat-1) call write_Igm_file(41,ispec2Dymax,NGLLX,NGLLZ,ilon,iz,izshift,ilayer_current)

           ! stzmin
           if (iz==0) then ! pas besoin du test comme précédemment car je stocke tout dans des tableaux et c'est pas
                             ! grave si on récrit les memes choses
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
               call write_Igm_file(42,ispec2Dzmin,NGLLX,NGLLY,ilon,ilat,0,ilayer_current)
              !open(125,file='ggl_elemts')
!!$              do kkk=1,1!NGLLZ
!!$                 do jjj=1,NGLLY
!!$                    do iii=1,NGLLX
!!$                       write(125,'(3f20.10)') xstore(iii,jjj,kkk)/1000.d0,  &
!!$                                       ystore(iii,jjj,kkk)/1000.d0,  zstore(iii,jjj,kkk)/1000.d0
!!$              !write(*,*) xstore
!!$                    enddo
!!$                 enddo
!!$              enddo
              !close(125)
              !   read(*,*) ia
              call store_zmin_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,&
                   lon_zmin,lat_zmin,nlon_dsm,nlat_dsm,ilon,ilat)
           endif

        enddo
     enddo
  enddo
  close(27)
  close(28)
  close(29)
  close(30)
  close(31)
  !close(32)

  ! ecriture des profondeurs de calcul pour DSM
  call write_recdepth_dsm(Ndepth,R_EARTH,MESH)
  ! ecriture de stzmin
  call write_stzmin(lon_zmin,lat_zmin,nlon_dsm,nlat_dsm,MESH)
  !

  z_bottom = minval(zgrid(:,:,:,:))
  zgrid(:,:,:,:) = zgrid(:,:,:,:) - z_bottom
  UTM_X_MIN=minval(xgrid)
  UTM_X_MAX=maxval(xgrid)
 ! modele 1D
  open(88,file=trim(MESH)//'model_1D.in')
  write(88,*) nlayer,4
  do i=1,nlayer
     write(88,*) zlayer(i)
     write(88,'(4f20.10)') vpv(i,:)
     write(88,'(4f20.10)') vsv(i,:)
     write(88,'(4f20.10)') density(i,:)
  enddo
  write(88,*)  z_bottom
  write(88,*)  lon_center_chunk,  lat_center_chunk,  chunk_azi
  close(88)
!
  !---------------- NUMEROTATION DES POINTS DE LA GRILLE ----


  ! on stocke touts les points de tous les elements
  do ispec=1,nspec

     ieoff = 8 * (ispec - 1)
     ilocnum = 0
     do k=1,2
        do j=1,2
           do i=1,2

              ilocnum = ilocnum + 1
              xp(ilocnum + ieoff)= xgrid(i,j,k,ispec)
              yp(ilocnum + ieoff)= ygrid(i,j,k,ispec)
              zp(ilocnum + ieoff)= zgrid(i,j,k,ispec)

           enddo
        enddo
     enddo
  enddo

  ! on identifie les points semblables et on les numerote
  call get_global1(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot,UTM_X_MIN,UTM_X_MAX)

  deallocate(xp,yp,zp)
  allocate(xp(nglob),yp(nglob),zp(nglob))

  ! on ne stocke que les points de la grille et leur numeros
  do ispec=1,nspec
     ieoff = 8 * (ispec - 1)
     ilocnum = 0
     do k=1,2
        do j=1,2
           do i=1,2
              ilocnum=ilocnum+1
              inum_loc(i,j,k,ispec) = iglob(ilocnum+ieoff)
              xp(iglob(ilocnum+ieoff)) = xgrid(i,j,k,ispec)
              yp(iglob(ilocnum+ieoff)) = ygrid(i,j,k,ispec)
              zp(iglob(ilocnum+ieoff)) = zgrid(i,j,k,ispec)
           enddo
        enddo
     enddo
  enddo

!---------------------------------------------------------------------

  write(90,*)  ispec2Dxmin
  write(90,*)  ispec2Dxmax
  write(90,*)  ispec2Dymin
  write(90,*)  ispec2Dymax
  write(90,*)  ispec2Dzmin
  close(27)
  close(28)
  close(29)
  close(30)
  close(31)
  close(32)
  close(37)
  close(38)
  close(39)
  close(40)
  close(41)
  close(42)
  close(81)
  close(82)
  close(83)
  close(84)
  close(85)
  close(86)
  close(87)
  close(88)
  close(89)
  close(90)
  !stop
  ! -------------------------------- SAUVEGARDE DES MESH FILES -----------

  open(27,file=trim(MESH)//'nodes_coords_file')
  write(27,*) nglob ! nb de sommets
  do kglob=1,nglob
     write(27,'(i14,3x,3(f20.5,1x))') kglob,xp(kglob),yp(kglob),zp(kglob)
  enddo
  close(27)

  open(27,file=trim(MESH)//'mesh_file')
  write(27,*) nspec
  do ispec=1,nspec
     write(27,'(9i15)')  ispec,inum_loc(1,1,1,ispec),inum_loc(2,1,1,ispec),&
          inum_loc(2,2,1,ispec),inum_loc(1,2,1,ispec),&
          inum_loc(1,1,2,ispec),inum_loc(2,1,2,ispec),&
          inum_loc(2,2,2,ispec),inum_loc(1,2,2,ispec)
  enddo
  close(27)
  !

  open(27,file=trim(MESH)//'absorbing_surface_file_xmin')
  write(27,*)  ispec2Dxmin
  do ispec=1,nspec
     if (iboun(1,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(1,1,1,ispec),inum_loc(1,2,1,ispec),&
          inum_loc(1,2,2,ispec),inum_loc(1,1,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_xmax')
  write(27,*) ispec2Dxmax
  do ispec=1,nspec
     if (iboun(2,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(2,1,1,ispec),inum_loc(2,2,1,ispec),&
          inum_loc(2,2,2,ispec),inum_loc(2,1,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_ymin')
  write(27,*) ispec2Dymin
  do ispec=1,nspec
     if (iboun(3,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(1,1,1,ispec),inum_loc(2,1,1,ispec),&
          inum_loc(2,1,2,ispec),inum_loc(1,1,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_ymax')
  write(27,*) ispec2Dymax
  do ispec=1,nspec
     if (iboun(4,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(1,2,1,ispec),inum_loc(2,2,1,ispec),&
          inum_loc(2,2,2,ispec),inum_loc(1,2,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_bottom')
  write(27,*) ispec2Dzmin
  do ispec=1,nspec
     if (iboun(5,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(1,1,1,ispec),inum_loc(1,2,1,ispec),&
          inum_loc(2,2,1,ispec),inum_loc(2,1,1,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'free_surface')
  write(27,*) ispec2Dzmax
  do ispec=1,nspec
     if (iboun(1,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(1,1,2,ispec),inum_loc(1,2,2,ispec),&
          inum_loc(2,2,2,ispec),inum_loc(2,1,2,ispec)
  enddo
  close(27)

  close(49)

  ! all processes done
  write(*,*) 'END '

  ! nothing to do anymore... bailing out
  stop

  end subroutine chunk_of_earth_Mesh

!=======================================================================================================
!
!=======================================================================================================

! compute the Euler angles and the associated rotation matrix

  subroutine euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

  implicit none

  !include "constants.h"

  double precision rotation_matrix(3,3)
  double precision CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH

  double precision alpha,beta,gamma
  double precision sina,cosa,sinb,cosb,sing,cosg

  double precision DEGREES_TO_RADIANS

  DEGREES_TO_RADIANS = 3.141592653589793d0/180.d0


! compute colatitude and longitude and convert to radians
  alpha = CENTER_LONGITUDE_IN_DEGREES * DEGREES_TO_RADIANS
  beta = (90.0d0 - CENTER_LATITUDE_IN_DEGREES) * DEGREES_TO_RADIANS
  gamma = GAMMA_ROTATION_AZIMUTH * DEGREES_TO_RADIANS

  sina = dsin(alpha)
  cosa = dcos(alpha)
  sinb = dsin(beta)
  cosb = dcos(beta)
  sing = dsin(gamma)
  cosg = dcos(gamma)

! define rotation matrix
  rotation_matrix(1,1) = cosg*cosb*cosa-sing*sina
  rotation_matrix(1,2) = -sing*cosb*cosa-cosg*sina
  rotation_matrix(1,3) = sinb*cosa
  rotation_matrix(2,1) = cosg*cosb*sina+sing*cosa
  rotation_matrix(2,2) = -sing*cosb*sina+cosg*cosa
  rotation_matrix(2,3) = sinb*sina
  rotation_matrix(3,1) = -cosg*sinb
  rotation_matrix(3,2) = sing*sinb
  rotation_matrix(3,3) = cosb

  end subroutine euler_angles

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_gllz_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,current_layer,nel_depth,ilayer,iz,Ndepth)

  implicit none

  integer NGLLX,NGLLY,NGLLZ,nel_depth,iz,Ndepth
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision profondeur
  integer current_layer(0:nel_depth-1),ilayer,k
  !write(*,*) ilayer,  current_layer(iz)
  !profondeur = dsqrt(xstore(1,1,k)**2 + ystore(1,1,k)**2 + (zstore(1,1,k) )**2 )
  !write(27,*) profondeur/1000., ilayer
  if (ilayer ==  current_layer(iz)) then
     do k=2,NGLLZ
        profondeur = dsqrt(xstore(1,1,k)**2 + ystore(1,1,k)**2 + (zstore(1,1,k) )**2 )
        write(27,*) profondeur/1000., ilayer-1,1
        Ndepth=Ndepth+1
     enddo
  else ! new layer

     k=1
     profondeur = dsqrt(xstore(1,1,k)**2 + ystore(1,1,k)**2 + (zstore(1,1,k) )**2 )
     if (ilayer==0) then
        ilayer =  current_layer(iz)
        write(27,*)  profondeur/1000., ilayer-1,1
        Ndepth=Ndepth+1
     else
        ilayer =  current_layer(iz)
        write(27,*)  profondeur/1000., ilayer-1,-1
        Ndepth=Ndepth+1
     endif
     do k=2,NGLLZ ! on duplique le dernier point
        profondeur = dsqrt(xstore(1,1,k)**2 + ystore(1,1,k)**2 + (zstore(1,1,k) )**2 )
        write(27,*)  profondeur/1000., ilayer-1,1
        Ndepth=Ndepth+1
     enddo


  endif

  end subroutine write_gllz_points

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_recdepth_dsm(Ndepth,R_EARTH,MESH)

  implicit none

  integer Ndepth,i
  double precision R_EARTH,prof
  double precision, allocatable :: z(:)
  integer, allocatable :: zindex(:),ziflag(:)
  integer ilayer,flag
  character(len=10) MESH



  open(27,file=trim(MESH)//'.recdepth')

  allocate(zindex(Ndepth),ziflag(Ndepth))
  allocate(z(Ndepth))

  do i=1,Ndepth
     read(27,*) prof,ilayer,flag
     z(Ndepth-i+1)=R_EARTH/1000.d0-prof
     zindex(Ndepth-i+1)=ilayer
     ziflag(Ndepth-i+1)=flag
  enddo

  close(27)
  open(27,file=trim(MESH)//'recdepth')

  write(27,*) Ndepth
  i=1
  write(27,*) z(i),zindex(i),ziflag(i)
  do i=2,Ndepth-1
     if (ziflag(i-1) == -1 ) then
        write(27,*) z(i),zindex(i),-1
     else
         write(27,*) z(i),zindex(i),1
     endif
  enddo
  i=Ndepth
  write(27,*) z(i),zindex(i),ziflag(i)

end subroutine write_recdepth_dsm

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stxmin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLY_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  if (test) then
     NGLLY_eff = NGLLY
  else
     NGLLY_eff = NGLLY - 1
  endif

  do jgll=1,NGLLY_eff
     vector_ori(1)=xstore(1,jgll,NGLLZ)
     vector_ori(2)=ystore(1,jgll,NGLLZ)
     vector_ori(3)=zstore(1,jgll,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      ! passage de geocentique à géographique
      !!theta = PI/2.D0 - lati
      ! convert the geocentric colatitude to a geographic colatitude
      !!colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
      !!lati =   PI/2.0d0 - colat

      !write(28,*) xstore(1,jgll,NGLLZ), ystore(1,jgll,NGLLZ), zstore(1,jgll,NGLLZ)!x,y !long/deg2rad,lati/deg2rad
      write(28,*) long/deg2rad,lati/deg2rad !,rayon/1000
      !write(38,'()') 1,(NGLLY-1)*jy_elm+jgll
       write(49,*)
       write(49,*)     vector_ori(:)
       write(49,*)     vector_rotated(:)

  enddo

  end subroutine write_stxmin

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stxmax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLY_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  if (test) then
     NGLLY_eff = NGLLY
  else
     NGLLY_eff = NGLLY - 1
  endif

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  do jgll=1,NGLLY_eff
     vector_ori(1)=xstore(NGLLX,jgll,NGLLZ)
     vector_ori(2)=ystore(NGLLX,jgll,NGLLZ)
     vector_ori(3) =zstore(NGLLX,jgll,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      ! passage de geocentique à géographique
      !!theta = PI/2.D0 - lati
      ! convert the geocentric colatitude to a geographic colatitude
      !!colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
      !!lati =   PI/2.0d0 - colat

      !write(28,*) xstore(1,jgll,NGLLZ), ystore(1,jgll,NGLLZ), zstore(1,jgll,NGLLZ)!x,y !long/deg2rad,lati/deg2rad
      write(29,*) long/deg2rad,lati/deg2rad !,rayon/1000
  enddo

  end subroutine write_stxmax

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stymin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLX_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

   if (test) then
     NGLLX_eff = NGLLX
  else
     NGLLX_eff = NGLLX - 1
  endif

  do jgll=1,NGLLX_eff
     vector_ori(1)=xstore(jgll,1,NGLLZ)
     vector_ori(2)=ystore(jgll,1,NGLLZ)
     vector_ori(3) =zstore(jgll,1,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      ! passage de geocentique à géographique
      !!theta = PI/2.D0 - lati
      ! convert the geocentric colatitude to a geographic colatitude
      !!colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
      !!lati =   PI/2.0d0 - colat

      !write(28,*) xstore(1,jgll,NGLLZ), ystore(1,jgll,NGLLZ), zstore(1,jgll,NGLLZ)!x,y !long/deg2rad,lati/deg2rad
      write(30,*) long/deg2rad,lati/deg2rad !,rayon/1000
  enddo

  end subroutine write_stymin

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stymax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLX_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  if (test) then
     NGLLX_eff = NGLLX
  else
     NGLLX_eff = NGLLX - 1
  endif

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  do jgll=1,NGLLX_eff
     vector_ori(1)=xstore(jgll,NGLLY,NGLLZ)
     vector_ori(2)=ystore(jgll,NGLLY,NGLLZ)
     vector_ori(3) =zstore(jgll,NGLLY,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      ! passage de geocentique à géographique
      !!theta = PI/2.D0 - lati
      ! convert the geocentric colatitude to a geographic colatitude
      !!colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
      !!lati =   PI/2.0d0 - colat

      !write(28,*) xstore(1,jgll,NGLLZ), ystore(1,jgll,NGLLZ), zstore(1,jgll,NGLLZ)!x,y !long/deg2rad,lati/deg2rad
      write(31,*) long/deg2rad,lati/deg2rad !,rayon/1000
  enddo

  end subroutine write_stymax

!=======================================================================================================
!
!=======================================================================================================

  subroutine store_zmin_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,&
             lon_zmin,lat_zmin,nlon_dsm,nlat_dsm,ilon,ilat)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,igll,jgll,i,j
  integer ilon,ilat,iglob,jglob,nlat_dsm,nlon_dsm
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  double precision lon_zmin(nlon_dsm,nlat_dsm),lat_zmin(nlon_dsm,nlat_dsm)


  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  do jgll=1,NGLLY
     do igll=1,NGLLX
        vector_ori(1)=xstore(igll,jgll,1)
        vector_ori(2)=ystore(igll,jgll,1)
        vector_ori(3) =zstore(igll,jgll,1)

        do i = 1,NDIM
           vector_rotated(i) = 0.d0
           do j = 1,NDIM
              vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
           enddo
        enddo
        x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
        rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

        long=atan2(y,x)
        lati=asin(z/rayon)

      ! passage de geocentique à géographique
      !!theta = PI/2.D0 - lati
      ! convert the geocentric colatitude to a geographic colatitude
      !!colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
      !!lati =   PI/2.0d0 - colat

      !write(28,*) xstore(1,jgll,NGLLZ), ystore(1,jgll,NGLLZ), zstore(1,jgll,NGLLZ)!x,y !long/deg2rad,lati/deg2rad
      !write(31,*) long/deg2rad,lati/deg2rad !,rayon/1000
        iglob=(ilon)*(NGLLX-1)+igll
        jglob=(ilat)*(NGLLY-1)+jgll
        lon_zmin(iglob,jglob)= long/deg2rad
        lat_zmin(iglob,jglob)= lati/deg2rad
        !write(32,'(3f20.10)') xstore(igll,jgll,1)/1000.d0, ystore(igll,jgll,1)/1000.d0,zstore(igll,jgll,1)/1000.d0
        !write(32,*) xstore(igll,jgll,NGLLZ), ystore(igll,igll,NGLLZ),zstore(igll,jgll,NGLLZ)
     enddo
  enddo

  end subroutine store_zmin_points

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stzmin(x,y,nx,ny,MESH)

  implicit none

  integer i,j,nx,ny
  double precision x(nx,ny),y(nx,ny)
  character(len=10) MESH

  open(27,file=trim(MESH)//'stzmin')
  write(27,*) nx*ny
  do j=1,ny
     do i=1,nx
        write(27,*) x(i,j),y(i,j)
     enddo
  enddo
  close(27)

  end subroutine write_stzmin

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_Igm_file(iunit,ispec2D,NGLL1,NGLL2,ie,je,js,il)

  implicit none

  integer iunit,ispec2D,NGLL1,NGLL2,ie,je,js,il
  integer i,j
  do j=1,NGLL2
     do i=1,NGLL1
        write(iunit,*) i,j,ispec2D,(NGLL1-1)*ie+i,(NGLL2-1)*je+j+js,il
     enddo
  enddo

  end subroutine write_Igm_file

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
!
!   ROUTINES POUR FAIRE DES ROTATIONS 3D ET DIVERS CHANGEMENTS DE REPERES
!
! Vadim Monteiller Mars 2013
!
!-------------------------------------------------------------------------------
! matrice de rotation 3D d'axe "axe" et d'angle theta (d°)
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

!------------------------------------------------------------------------------
! rotation pour passer d'un repere local a un autre
