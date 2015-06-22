!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
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

  subroutine earth_chunk_HEX27_Mesher(NGNOD)

  use constants, only: NGLLX, NGLLY, NGLLZ, NDIM, R_EARTH, PI, ZERO, TINYVAL, &
                       old_DSM_coupling_from_Vadim, EXTERNAL_CODE_IS_AXISEM, EXTERNAL_CODE_IS_DSM

  use shared_parameters, only: EXTERNAL_CODE_TYPE

  implicit none

!==============================================================================================!
!                                                                                              !
!  Singular option of meshfem3D : MESH OF A GLOBE EARTH CHUNK FOR THE INTERFACE DSM-SPECFEM3D  !
!  Case of 27 nodes per element (HEX27)                                                        !
!                                                                                              !
!  Integrated in meshfem3d by CD, October 2014                                                 !
!                                                                                              !
!  WARNING : A local convention is used for the mapping of                                     !
!            the cubic sphere (to complete)                                                    !
!                                                                                              !
!==============================================================================================!

!
!--- Parameters
!

  integer, parameter :: myrank = 0
  integer, parameter :: nlayer = 12 !! (number of layer in the model iasp91, or ak135, or prem (one more layer than the model)

  double precision, parameter :: GAUSSALPHA = 0.d0, GAUSSBETA = 0.d0

  logical, parameter ::  RUN_BENCHMARK = .false.

!
!--- Other
!

  integer NGNOD

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
  double precision deg2rad
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

  integer ::  istore_for_new_outputs
  integer ::   updown(NGLLZ)
  double precision , dimension(NGLLX,NGLLY,NGLLZ) ::  longitud, latitud, radius

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
!---                          (k = 6 with -z for the mapping of the cubic sphere, cf Chevrot 2012)
!---                          We define the mesh of a chunk of the earth in the cubic sphere
!

  deg2rad = 3.141592653589793d0/180.d0

  open(49, file=trim(MESH)//'output_mesher_chunk_HEX27.txt')

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
  nglob    = (2*nel_lat + 1) * (2*nel_lon + 1) * (2*nel_depth + 1)
  nspec    = nel_lat * nel_lon * nel_depth
  npointot = 27 * nspec

  allocate(xp(npointot), yp(npointot), zp(npointot))
  allocate(iglob(npointot), loc(npointot))
  allocate(ifseg(npointot))
  allocate(ProfForGemini(0:NZ-1,3))
  allocate(current_layer(0:NZ-1))
  allocate(inum_loc(3,3,3,nspec))
  allocate(xgrid(3,3,3,nspec), ygrid(3,3,3,nspec), zgrid(3,3,3,nspec))
  allocate(lon_zmin(nlon_dsm,nlat_dsm), lat_zmin(nlon_dsm,nlat_dsm))
  allocate(iboun(6,nspec)) ! boundary locator

  iboun(:,:) = .false.


!! MODIF HEX27 LA ==> cf call hex_nodes-----------------------------

 ! corner nodes

  iaddx(1) = 0
  iaddy(1) = 0
  iaddz(1) = 0

  iaddx(2) = 2
  iaddy(2) = 0
  iaddz(2) = 0

  iaddx(3) = 2
  iaddy(3) = 2
  iaddz(3) = 0

  iaddx(4) = 0
  iaddy(4) = 2
  iaddz(4) = 0

  iaddx(5) = 0
  iaddy(5) = 0
  iaddz(5) = 2

  iaddx(6) = 2
  iaddy(6) = 0
  iaddz(6) = 2

  iaddx(7) = 2
  iaddy(7) = 2
  iaddz(7) = 2

  iaddx(8) = 0
  iaddy(8) = 2
  iaddz(8) = 2

! midside nodes (nodes located in the middle of an edge)

  iaddx(9) = 1
  iaddy(9) = 0
  iaddz(9) = 0

  iaddx(10) = 2
  iaddy(10) = 1
  iaddz(10) = 0

  iaddx(11) = 1
  iaddy(11) = 2
  iaddz(11) = 0

  iaddx(12) = 0
  iaddy(12) = 1
  iaddz(12) = 0

  iaddx(13) = 0
  iaddy(13) = 0
  iaddz(13) = 1

  iaddx(14) = 2
  iaddy(14) = 0
  iaddz(14) = 1

  iaddx(15) = 2
  iaddy(15) = 2
  iaddz(15) = 1

  iaddx(16) = 0
  iaddy(16) = 2
  iaddz(16) = 1

  iaddx(17) = 1
  iaddy(17) = 0
  iaddz(17) = 2

  iaddx(18) = 2
  iaddy(18) = 1
  iaddz(18) = 2

  iaddx(19) = 1
  iaddy(19) = 2
  iaddz(19) = 2

  iaddx(20) = 0
  iaddy(20) = 1
  iaddz(20) = 2

! side center nodes (nodes located in the middle of a face)

  iaddx(21) = 1
  iaddy(21) = 1
  iaddz(21) = 0

  iaddx(22) = 1
  iaddy(22) = 0
  iaddz(22) = 1

  iaddx(23) = 2
  iaddy(23) = 1
  iaddz(23) = 1

  iaddx(24) = 1
  iaddy(24) = 2
  iaddz(24) = 1

  iaddx(25) = 0
  iaddy(25) = 1
  iaddz(25) = 1

  iaddx(26) = 1
  iaddy(26) = 1
  iaddz(26) = 2

! center node (barycenter of the eight corners)

  iaddx(27) = 1
  iaddy(27) = 1
  iaddz(27) = 1

!! --------------------------------------------
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

!
!--- for new output mesh files and VM coupling with AxiSEM
!

  open(91, file = trim(MESH)//'list_ggl_boundary_spherical.txt')
  open(92, file = trim(MESH)//'list_ggl_boundary_cartesian.txt')

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

           istore_for_new_outputs = 0

           ! get boundary

           ! on boundary 1: x=xmin
           if (ilon == 0) then

              iboun(1,ispec)=.true.
              ispec2Dxmin=ispec2Dxmin+1
              write(89,*) ispec,ispec2Dxmin,1

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 2: xmax
           if (ilon == nel_lon-1) then

              iboun(2,ispec)=.true.
              ispec2Dxmax=ispec2Dxmax+1
              !write(*,*) '------ TOZ',ispec,ilon
              write(89,*) ispec,ispec2Dxmax,2

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 3: ymin
           if (ilat == 0) then

              iboun(3,ispec)=.true.
              ispec2Dymin=ispec2Dymin+1
              write(89,*) ispec,ispec2Dymin,3

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 4: ymax
           if (ilat == nel_lat-1) then

              iboun(4,ispec) =.true.
              ispec2Dymax=ispec2Dymax+1
              write(89,*) ispec,ispec2Dymax,4

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 5: bottom
           if (iz == 0) then

              iboun(5,ispec)=.true.
              ispec2Dzmin=ispec2Dzmin+1
              write(89,*) ispec,ispec2Dzmin,5

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 6: top
           if (iz == nel_depth-1) then
              ispec2Dzmax= ispec2Dzmax+1
              iboun(6,ispec)=.true.
           endif

           do ia=1,NGNOD

!! MODIF HEX27 LA -----------------------------

              i=iaddx(ia)
              j=iaddy(ia)
              k=iaddz(ia)

              SELECT CASE (k)
                CASE(0)
                  z = 1000d0*ProfForGemini(iz,1)
                CASE(1)
                  z = 1000d0*ProfForGemini(iz,3)
                CASE(2)
                  z = 1000d0*ProfForGemini(iz,2)
              END SELECT

              ! longitude
              ratio_xi = (dble(ilon) + dble(i)/2.d0 ) / dble(NX)
              x = 2.d0*ratio_xi-1.d0
              x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

              ! latitude
              ratio_eta = (dble(ilat) + dble(j)/2.d0 ) / dble(NY)
              y = 2.d0*ratio_eta-1.d0
              y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

              ! mapping cubic sphere (k=6, Chevrot at al 2012, avec -z)

              pz= z/dsqrt(1.d0 + y*y + x*x) !(=r/s)
              px= pz * x !(tan(xi) * r/s)
              py= pz * y !(tan(eta) * r/s)

              ! old version
              xgrid(i+1,j+1,k+1,ispec) = px !px
              ygrid(i+1,j+1,k+1,ispec) = py !py
              zgrid(i+1,j+1,k+1,ispec) = pz

              xelm(ia)=xgrid(i+1,j+1,k+1,ispec)
              yelm(ia)=ygrid(i+1,j+1,k+1,ispec)
              zelm(ia)=zgrid(i+1,j+1,k+1,ispec)

           enddo

           ! INTERFACE FOR DSM ------

           ! Vertical receptors

           if (ilat==0 .and. ilon==0) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              call write_gllz_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,current_layer,nel_depth,ilayer,iz,Ndepth,updown)
           endif

          ! Write two files giving Spherical coordinate on ALL the GLL points on the surface of the 3D chunk for the new DSM
          ! coupling (light version using 2D chunk)
          !
          ! ==> CAUTION : will be also used later for the VM coupling with AxiSEM
          !
          ! (must be after write_gllz_points to know the value of ilayer)

          if ( ( ( EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_DSM .and. (.not. old_DSM_coupling_from_Vadim) ) .or.    &
                 ( EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM                                        ) ) .and. &
               ( istore_for_new_outputs > 0 ) ) then


            call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
            call write_all_chunk_surface_GLL_in_spherical_and_cartesian_coords(xstore,ystore,zstore, &
                                                          deg2rad,ilayer,iboun,ispec,nspec,longitud, &
                                                          latitud,radius,rotation_matrix,updown)

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
           if (iz==0) then ! pas besoin du test comme precedemment car je stocke tout dans des tableaux et c'est pas
                             ! grave si on recrit les memes choses
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)

              call write_Igm_file(42,ispec2Dzmin,NGLLX,NGLLY,ilon,ilat,0,ilayer_current)

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

     ieoff = 27 * (ispec - 1)
     ilocnum = 0

     do k=1,3
        do j=1,3
           do i=1,3

              ilocnum = ilocnum + 1
              xp(ilocnum + ieoff)= xgrid(i,j,k,ispec)
              yp(ilocnum + ieoff)= ygrid(i,j,k,ispec)
              zp(ilocnum + ieoff)= zgrid(i,j,k,ispec)

           enddo
        enddo
     enddo
  enddo

  ! on identifie les points semblables et on les numerote
  call get_global1(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot,NGNOD,UTM_X_MIN,UTM_X_MAX)

  deallocate(xp,yp,zp)
  allocate(xp(nglob),yp(nglob),zp(nglob))

!! MODIF HEX27 LA -----------------------------

  ! on ne stocke que les points de la grille et leur numeros
  do ispec=1,nspec

     ieoff = 27 * (ispec - 1)
     ilocnum = 0

     do k=1,3
        do j=1,3
           do i=1,3

              ilocnum                  = ilocnum + 1
              inum_loc(i,j,k,ispec)    = iglob(ilocnum+ieoff)
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
  close(91)
  close(92)
  !stop
  ! -------------------------------- SAUVEGARDE DES MESH FILES -----------

!! MODIF HEX27 LA -----------------------------

  open(27,file=trim(MESH)//'nodes_coords_file')
  write(27,*) nglob ! nb de sommets
  do kglob=1,nglob
     write(27,'(i14,3x,3(f20.5,1x))') kglob,xp(kglob),yp(kglob),zp(kglob)
  enddo
  close(27)

  open(27,file=trim(MESH)//'mesh_file')
  write(27,*) nspec
  do ispec=1,nspec
    write(27,'(28i15)') ispec, &
                        inum_loc(1,1,1,ispec), inum_loc(3,1,1,ispec), inum_loc(3,3,1,ispec), inum_loc(1,3,1,ispec), &
                        inum_loc(1,1,3,ispec), inum_loc(3,1,3,ispec), inum_loc(3,3,3,ispec), inum_loc(1,3,3,ispec), &
                        inum_loc(2,1,1,ispec), inum_loc(3,2,1,ispec), inum_loc(2,3,1,ispec), inum_loc(1,2,1,ispec), &
                        inum_loc(1,1,2,ispec), inum_loc(3,1,2,ispec), inum_loc(3,3,2,ispec), inum_loc(1,3,2,ispec), &
                        inum_loc(2,1,3,ispec), inum_loc(3,2,3,ispec), inum_loc(2,3,3,ispec), inum_loc(1,2,3,ispec), &
                        inum_loc(2,2,1,ispec), inum_loc(2,1,2,ispec), inum_loc(3,2,2,ispec), inum_loc(2,3,2,ispec), &
                        inum_loc(1,2,2,ispec), inum_loc(2,2,3,ispec), inum_loc(2,2,2,ispec)
  enddo
  close(27)
  !

  open(27,file=trim(MESH)//'absorbing_surface_file_xmin')
  write(27,*)  ispec2Dxmin
  do ispec=1,nspec
     if (iboun(1,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(1,1,1,ispec), inum_loc(1,3,1,ispec), &
                                                  inum_loc(1,3,3,ispec), inum_loc(1,1,3,ispec), &
                                                  inum_loc(1,2,1,ispec), inum_loc(1,3,2,ispec), &
                                                  inum_loc(1,2,3,ispec), inum_loc(1,1,2,ispec), &
                                                  inum_loc(1,2,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_xmax')
  write(27,*) ispec2Dxmax
  do ispec=1,nspec
     if (iboun(2,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(3,1,1,ispec), inum_loc(3,3,1,ispec), &
                                                  inum_loc(3,3,3,ispec), inum_loc(3,1,3,ispec), &
                                                  inum_loc(3,2,1,ispec), inum_loc(3,3,2,ispec), &
                                                  inum_loc(3,2,3,ispec), inum_loc(3,1,2,ispec), &
                                                  inum_loc(3,2,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_ymin')
  write(27,*) ispec2Dymin
  do ispec=1,nspec
     if (iboun(3,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(1,1,1,ispec), inum_loc(3,1,1,ispec), &
                                                  inum_loc(3,1,3,ispec), inum_loc(1,1,3,ispec), &
                                                  inum_loc(2,1,1,ispec), inum_loc(3,1,2,ispec), &
                                                  inum_loc(2,1,3,ispec), inum_loc(1,1,2,ispec), &
                                                  inum_loc(2,1,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_ymax')
  write(27,*) ispec2Dymax
  do ispec=1,nspec
     if (iboun(4,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(1,3,1,ispec), inum_loc(3,3,1,ispec), &
                                                  inum_loc(3,3,3,ispec), inum_loc(1,3,3,ispec), &
                                                  inum_loc(2,3,1,ispec), inum_loc(3,3,2,ispec), &
                                                  inum_loc(2,3,3,ispec), inum_loc(1,3,2,ispec), &
                                                  inum_loc(2,3,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_bottom')
  write(27,*) ispec2Dzmin
  do ispec=1,nspec
     if (iboun(5,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(1,1,1,ispec), inum_loc(3,1,1,ispec), &
                                                  inum_loc(3,3,1,ispec), inum_loc(1,3,1,ispec), &
                                                  inum_loc(2,1,1,ispec), inum_loc(3,2,1,ispec), &
                                                  inum_loc(2,3,1,ispec), inum_loc(1,2,1,ispec), &
                                                  inum_loc(2,2,1,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'free_surface')
  write(27,*) ispec2Dzmax
  do ispec=1,nspec
     if (iboun(6,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(1,1,3,ispec), inum_loc(3,1,3,ispec), &
                                                  inum_loc(3,3,3,ispec), inum_loc(1,3,3,ispec), &
                                                  inum_loc(2,1,3,ispec), inum_loc(3,2,3,ispec), &
                                                  inum_loc(2,3,3,ispec), inum_loc(1,2,3,ispec), &
                                                  inum_loc(2,2,3,ispec)
  enddo
  close(27)

  close(49)

  ! all processes done
  write(*,*) 'END '

  end subroutine earth_chunk_HEX27_Mesher

!=======================================================================================================!
