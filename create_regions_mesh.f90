!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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


  subroutine create_regions_mesh_ext_mesh(ibool, &
           xstore,ystore,zstore,nspec,npointot,myrank,LOCAL_PATH, &
           nnodes_ext_mesh,nelmnts_ext_mesh, &
           nodes_coords_ext_mesh,elmnts_ext_mesh,max_static_memory_size,mat_ext_mesh,materials_ext_mesh, &
           nmat_ext_mesh,undef_mat_prop,nundefMat_ext_mesh,ninterface_ext_mesh,max_interface_size_ext_mesh, &
           my_neighbours_ext_mesh,my_nelmnts_neighbours_ext_mesh,my_interfaces_ext_mesh, &
           ibool_interfaces_ext_mesh,nibool_interfaces_ext_mesh, &
           nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, NSPEC2D_BOTTOM, NSPEC2D_TOP,&
           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
           ibelm_xmin,ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top)

! create the different regions of the mesh

  implicit none

  include "constants.h"

! number of spectral elements in each block
  integer nspec

  integer npointot

  character(len=150) LOCAL_PATH

! arrays with the mesh
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

! static memory size needed by the solver
  double precision :: static_memory_size,max_static_memory_size 

! data from the external mesh
  integer :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  integer, dimension(ESIZE,nelmnts_ext_mesh) :: elmnts_ext_mesh
  integer, dimension(2,nelmnts_ext_mesh) :: mat_ext_mesh
!  double precision, external :: materials_ext_mesh
  integer :: ninterface_ext_mesh,max_interface_size_ext_mesh
  integer, dimension(ninterface_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(ninterface_ext_mesh) :: my_nelmnts_neighbours_ext_mesh
  integer, dimension(6,max_interface_size_ext_mesh,ninterface_ext_mesh) :: my_interfaces_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,ninterface_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(ninterface_ext_mesh) :: nibool_interfaces_ext_mesh
!pll
  integer :: nmat_ext_mesh,nundefMat_ext_mesh 
  double precision, dimension(5,nmat_ext_mesh) :: materials_ext_mesh  
  character (len=30), dimension(5,nundefMat_ext_mesh):: undef_mat_prop

! for MPI buffers
  integer, dimension(:), allocatable :: reorder_interface_ext_mesh,ind_ext_mesh,ninseg_ext_mesh,iwork_ext_mesh
  integer, dimension(:), allocatable :: nibool_interfaces_ext_mesh_true
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer, dimension(:), allocatable :: ibool_interface_ext_mesh_dummy
  double precision, dimension(:), allocatable :: work_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore_dummy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: ystore_dummy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: zstore_dummy

! Gauss-Lobatto-Legendre points and weights of integration
  double precision, dimension(:), allocatable :: xigll,yigll,zigll,wxgll,wygll,wzgll

! 3D shape functions and their derivatives
  double precision, dimension(:,:,:,:), allocatable :: shape3D
  double precision, dimension(:,:,:,:,:), allocatable :: dershape3D

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

! the jacobian
  real(kind=CUSTOM_REAL) jacobianl

! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore

! for model density
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore,kappastore,mustore,vpstore,vsstore 

! attenuation 
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: iflag_attenuation_store

! proc numbers for MPI
  integer myrank

! check area and volume of the final mesh
  double precision weight

! variables for creating array ibool (some arrays also used for AVS or DX files)
  integer, dimension(:), allocatable :: iglob,locval
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,zp

  integer nglob
  integer ieoff,ilocnum
  integer ier
  integer iinterface

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass

! ---------------------------

! name of the database file
  character(len=150) prname

  integer i,j,k,ia,ispec,iglobnum

! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  integer, dimension(:,:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber

! pll 
  integer :: iundef
  logical, dimension(6,nspec) :: iboun  

! 2D shape functions and their derivatives
  double precision, dimension(:,:,:), allocatable :: shape2D_x,shape2D_y,shape2D_bottom,shape2D_top
  double precision, dimension(:,:,:,:), allocatable :: dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top

  integer  :: ispec2D
  integer :: iflag, flag_below, flag_above
  integer  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, NSPEC2D_BOTTOM, NSPEC2D_TOP
  integer  :: NSPEC2DMAX_XMIN_XMAX 
  integer  :: NSPEC2DMAX_YMIN_YMAX
  integer, dimension(nspec2D_xmin)  :: ibelm_xmin  
  integer, dimension(nspec2D_xmax)  :: ibelm_xmax
  integer, dimension(nspec2D_ymin)  :: ibelm_ymin
  integer, dimension(nspec2D_ymax)  :: ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM)  :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP)  :: ibelm_top
  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX) :: nimin,nimax
  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX) :: njmin,njmax 
  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX) :: nkmin_xi
  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX) :: nkmin_eta

  ! 2-D jacobians and normals
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
       jacobian2D_xmin,jacobian2D_xmax, &
       jacobian2D_ymin,jacobian2D_ymax,jacobian2D_bottom,jacobian2D_top
  
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top
  
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs



! For Piero Basini :
! integer :: doubling_value_found_for_Piero
!   double precision :: xmesh,ymesh,zmesh
!   double precision :: rho,vp,vs

!   integer,dimension(nspec) ::  idoubling
!   integer :: doubling_value_found_for_Piero
!   integer, parameter :: NUMBER_OF_STATIONS = 6
!   double precision, parameter :: RADIUS_TO_EXCLUDE = 250.d0
!   double precision, dimension(NUMBER_OF_STATIONS) :: utm_x_station,utm_y_station

!   logical :: is_around_a_station
!   integer :: istation

! ! store bedrock values
!   integer ::  icornerlat,icornerlong
!   double precision ::  lat,long,elevation_bedrock
!   double precision ::  lat_corner,long_corner,ratio_xi,ratio_eta
  
! ! size of topography and bathymetry file for Piero Basini's model
!   integer, parameter :: NX_TOPO = 787, NY_TOPO = 793
!   double precision, parameter :: ORIG_LAT_TOPO = -102352.d0
!   double precision, parameter :: ORIG_LONG_TOPO = 729806.d0
!   character(len=150), parameter :: TOPO_FILE = 'DATA/piero_model/dem_EV_UTM_regular_250_reordered.dat'
! ! for Piero Basini's model this is the resolution in meters of the topo file
!   double precision, parameter :: DEGREES_PER_CELL_TOPO = 250.d0

!real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: ibedrock


! **************

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

! Gauss-Lobatto-Legendre points of integration
  allocate(xigll(NGLLX))
  allocate(yigll(NGLLY))
  allocate(zigll(NGLLZ))

! Gauss-Lobatto-Legendre weights of integration
  allocate(wxgll(NGLLX))
  allocate(wygll(NGLLY))
  allocate(wzgll(NGLLZ))

! 3D shape functions and their derivatives
  allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ))
  allocate(dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ))

! pll 2D shape functions and their derivatives
  allocate(shape2D_x(NGNOD2D,NGLLY,NGLLZ))
  allocate(shape2D_y(NGNOD2D,NGLLX,NGLLZ))
  allocate(shape2D_bottom(NGNOD2D,NGLLX,NGLLY))
  allocate(shape2D_top(NGNOD2D,NGLLX,NGLLY))
  allocate(dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ))
  allocate(dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ))
  allocate(dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY))
  allocate(dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY))

! pll Stacey
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho_vs(NGLLX,NGLLY,NGLLZ,nspec))

! array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(mustore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vpstore(NGLLX,NGLLY,NGLLZ,nspec)) !pll
  allocate(vsstore(NGLLX,NGLLY,NGLLZ,nspec)) !pll

! arrays with mesh parameters
  allocate(xixstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xiystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xizstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etaxstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etaystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etazstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammaxstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammaystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammazstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(jacobianstore(NGLLX,NGLLY,NGLLZ,nspec))

! pll 2-D jacobians and normals
  allocate(jacobian2D_xmin(NGLLY,NGLLZ,nspec2D_xmin))
  allocate(jacobian2D_xmax(NGLLY,NGLLZ,nspec2D_xmax))
  allocate(jacobian2D_ymin(NGLLX,NGLLZ,nspec2D_ymin))
  allocate(jacobian2D_ymax(NGLLX,NGLLZ,nspec2D_ymax))
  allocate(jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM))
  allocate(jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP))

  allocate(normal_xmin(NDIM,NGLLY,NGLLZ,nspec2D_xmin))
  allocate(normal_xmax(NDIM,NGLLY,NGLLZ,nspec2D_xmax))
  allocate(normal_ymin(NDIM,NGLLX,NGLLZ,nspec2D_ymin))
  allocate(normal_ymax(NDIM,NGLLX,NGLLZ,nspec2D_ymax))
  allocate(normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM))
  allocate(normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP))

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly zero
  if(mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = ZERO
  if(mod(NGLLY,2) /= 0) yigll((NGLLY-1)/2+1) = ZERO
  if(mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = ZERO

! get the 3-D shape functions
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll)

! pll get the 2-D shape functions
  call get_shape2D(myrank,shape2D_x,dershape2D_x,yigll,zigll,NGLLY,NGLLZ)
  call get_shape2D(myrank,shape2D_y,dershape2D_y,xigll,zigll,NGLLX,NGLLZ)
  call get_shape2D(myrank,shape2D_bottom,dershape2D_bottom,xigll,yigll,NGLLX,NGLLY)
  call get_shape2D(myrank,shape2D_top,dershape2D_top,xigll,yigll,NGLLX,NGLLY)

! allocate memory for arrays
  allocate(iglob(npointot))
  allocate(locval(npointot))
  allocate(ifseg(npointot))
  allocate(xp(npointot))
  allocate(yp(npointot))
  allocate(zp(npointot))

!---

  xstore(:,:,:,:) = 0.d0
  ystore(:,:,:,:) = 0.d0
  zstore(:,:,:,:) = 0.d0

  do ispec = 1, nspec
     !call get_xyzelm(xelm, yelm, zelm, ispec, elmnts_ext_mesh, nodes_coords_ext_mesh, nspec, nnodes_ext_mesh)
     do ia = 1,NGNOD
     xelm(ia) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(ia,ispec))
     yelm(ia) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(ia,ispec))
     zelm(ia) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(ia,ispec))
     enddo

     call calc_jacobian(myrank,xixstore,xiystore,xizstore, &
          etaxstore,etaystore,etazstore, &
          gammaxstore,gammaystore,gammazstore,jacobianstore, &
          xstore,ystore,zstore, &
          xelm,yelm,zelm,shape3D,dershape3D,ispec,nspec)

  enddo

! !  Piero, read bedrock file
!  allocate(ibedrock(NX_TOPO_ANT,NY_TOPO_ANT))              
!  if(myrank == 0) then
!      call read_bedrock_file(ibedrock)
!  !    write(IMAIN,*)
!  !    write(IMAIN,*) 'regional bedrock file read ranges in m from ',minval(ibedrock),' to ',maxval(ibedrock)
!  !    write(IMAIN,*)
!   endif
!  ! broadcast the information read on the master to the nodes
!  ! call MPI_BCAST(ibedrock,NX_TOPO_ANT*NY_TOPO_ANT,MPI_REAL,0,MPI_COMM_WORLD,ier)
! call bcast_all_cr(ibedrock,NX_TOPO_ANT*NY_TOPO_ANT)

! kappastore and mustore
  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
! check if the material is known or unknown
           if (mat_ext_mesh(1,ispec) > 0) then
              rhostore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(1,ispec))
              vpstore(i,j,k,ispec) = materials_ext_mesh(2,mat_ext_mesh(1,ispec))
              vsstore(i,j,k,ispec) = materials_ext_mesh(3,mat_ext_mesh(1,ispec))
              iflag_attenuation_store(i,j,k,ispec) = materials_ext_mesh(4,mat_ext_mesh(1,ispec))
              !change for piero :
              !if(mat_ext_mesh(1,ispec) == 1) then
              !   iflag_attenuation_store(i,j,k,ispec) = 1
              !else
              !   iflag_attenuation_store(i,j,k,ispec) = 2
              !endif
           else if (mat_ext_mesh(2,ispec) == 1) then
              do iundef = 1,nundefMat_ext_mesh 
                 if(trim(undef_mat_prop(2,iundef)) == 'interface') then
                    read(undef_mat_prop(4,iundef),'(1i3)') flag_below
                    read(undef_mat_prop(5,iundef),'(1i3)') flag_above
                 endif
              enddo
              !call interface(iflag,flag_below,flag_above,ispec,nspec,i,j,k,xstore,ystore,zstore,ibedrock)
              rhostore(i,j,k,ispec) = materials_ext_mesh(1,iflag)
              vpstore(i,j,k,ispec) = materials_ext_mesh(2,iflag)
              vsstore(i,j,k,ispec) = materials_ext_mesh(3,iflag)
              iflag_attenuation_store(i,j,k,ispec) = materials_ext_mesh(4,iflag)
              !change for piero :
              !  if(iflag == 1) then
              !     iflag_attenuation_store(i,j,k,ispec) = 1
              !  else
              !     iflag_attenuation_store(i,j,k,ispec) = 2
              !  endif
             else
             ! call tomography()
           end if

           kappastore(i,j,k,ispec) = rhostore(i,j,k,ispec)*(vpstore(i,j,k,ispec)*vpstore(i,j,k,ispec) - &
                4.d0*vsstore(i,j,k,ispec)*vsstore(i,j,k,ispec)/3.d0)
           mustore(i,j,k,ispec) = rhostore(i,j,k,ispec)*vsstore(i,j,k,ispec)*&
                vsstore(i,j,k,ispec)
           
           ! Stacey, a completer par la suite  
           rho_vp(i,j,k,ispec) = rhostore(i,j,k,ispec)*vpstore(i,j,k,ispec)
           rho_vs(i,j,k,ispec) = rhostore(i,j,k,ispec)*vsstore(i,j,k,ispec)
           !end pll

        enddo
      enddo
    enddo
  enddo


! !! DK DK store the position of the six stations to be able to
! !! DK DK exclude circles around each station to make sure they are on the bedrock
! !! DK DK and not in the ice
!    utm_x_station(1) =  783500.6250000d0
!    utm_y_station(1) = -11828.7519531d0

!    utm_x_station(2) =  853644.5000000d0
!    utm_y_station(2) = -114.0138092d0

!    utm_x_station(3) = 863406.0000000d0
!    utm_y_station(3) = -53736.1640625d0

!    utm_x_station(4) =   823398.8125000d0
!    utm_y_station(4) = 29847.4511719d0

!    utm_x_station(5) = 863545.3750000d0
!    utm_y_station(5) = 19669.6621094d0

!    utm_x_station(6) =  817099.3750000d0
!    utm_y_station(6) = -24430.2871094d0

!  print*,myrank,'après store the position of the six stations'
!  call flush(6)

!  print*, myrank,minval(nodes_coords_ext_mesh(1,:))
!  call flush(6)


! print*, myrank,maxval(nodes_coords_ext_mesh(1,:))
!  call flush(6)


!  do ispec = 1, nspec

!     zmesh = zstore(2,2,2,ispec)

!    ! if(doubling_index == IFLAG_ONE_LAYER_TOPOGRAPHY) then
!     if(any(ibelm_top == ispec)) then
!     doubling_value_found_for_Piero = IFLAG_ONE_LAYER_TOPOGRAPHY
       
!     else if(zmesh < Z_23p4km) then
!        doubling_value_found_for_Piero = IFLAG_MANTLE_BELOW_23p4km
       
!     else if(zmesh < Z_14km) then
!        doubling_value_found_for_Piero = IFLAG_14km_to_23p4km
       
!     else
!        doubling_value_found_for_Piero = IFLAG_BEDROCK_down_to_14km
!     endif
!    idoubling(ispec) = doubling_value_found_for_Piero

!     do k = 1, NGLLZ
!       do j = 1, NGLLY
!         do i = 1, NGLLX

           
!            if(idoubling(ispec) == IFLAG_ONE_LAYER_TOPOGRAPHY .or. idoubling(ispec) == IFLAG_BEDROCK_down_to_14km) then
              
!               ! since we have suppressed UTM projection for Piero Basini, UTMx is the same as long
!               ! and UTMy is the same as lat
!               long = xstore(i,j,k,ispec)
!               lat = ystore(i,j,k,ispec)
              
!               ! get coordinate of corner in model
!               icornerlong = int((long - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
!               icornerlat = int((lat - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1
              
!               ! avoid edge effects and extend with identical point if outside model
!               if(icornerlong < 1) icornerlong = 1
!               if(icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
!               if(icornerlat < 1) icornerlat = 1
!               if(icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1
              
!               ! compute coordinates of corner
!               long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
!               lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO
                   
!               ! compute ratio for interpolation
!               ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO
!               ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO
                   
!               ! avoid edge effects
!               if(ratio_xi < 0.) ratio_xi = 0.
!               if(ratio_xi > 1.) ratio_xi = 1.
!               if(ratio_eta < 0.) ratio_eta = 0.
!               if(ratio_eta > 1.) ratio_eta = 1.
                   
!               ! interpolate elevation at current point
!               elevation_bedrock = &
!                    ibedrock(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
!                    ibedrock(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
!                    ibedrock(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
!                    ibedrock(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta
                   
!               !! DK DK exclude circles around each station to make sure they are on the bedrock
!               !! DK DK and not in the ice
!               is_around_a_station = .false.
!               do istation = 1,NUMBER_OF_STATIONS
!                  if(sqrt((long - utm_x_station(istation))**2 + (lat - utm_y_station(istation))**2) < RADIUS_TO_EXCLUDE) then
!                     is_around_a_station = .true.
!                     exit
!                  endif
!               enddo
              
!               ! define elastic parameters in the model
              
!               ! we are above the bedrock interface i.e. in the ice, and not too close to a station
!               if(zmesh >= elevation_bedrock .and. .not. is_around_a_station) then
!                  vp = 3800.d0
!                  vs = 1900.d0
!                  rho = 900.d0
!                  iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_ICE
                 
!                  ! we are below the bedrock interface i.e. in the bedrock, or close to a station
!               else
!                  vp = 5800.d0
!                  vs = 3200.d0
!                  rho = 2600.d0
!                  iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_BEDROCK
!               endif
              
!            else if(idoubling(ispec) == IFLAG_14km_to_23p4km) then
!               vp = 6800.d0
!               vs = 3900.d0
!               rho = 2900.d0
!               iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_BEDROCK
              
!            else if(idoubling(ispec) == IFLAG_MANTLE_BELOW_23p4km) then
!               vp = 8100.d0
!               vs = 4480.d0
!               rho = 3380.d0
!               iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_BEDROCK
              
!            endif
           
!                 !pll  8/06
!                     if(CUSTOM_REAL == SIZE_REAL) then
!                        rhostore(i,j,k,ispec) = sngl(rho)
!                        vpstore(i,j,k,ispec) = sngl(vp)
!                        vsstore(i,j,k,ispec) = sngl(vs)
!                     else
!                        rhostore(i,j,k,ispec) = rho
!                        vpstore(i,j,k,ispec) = vp
!                        vsstore(i,j,k,ispec) = vs
!                     end if
                
!                 kappastore(i,j,k,ispec) = rhostore(i,j,k,ispec)*(vpstore(i,j,k,ispec)*vpstore(i,j,k,ispec) - &
!                      4.d0*vsstore(i,j,k,ispec)*vsstore(i,j,k,ispec)/3.d0)
!                 mustore(i,j,k,ispec) = rhostore(i,j,k,ispec)*vsstore(i,j,k,ispec)*&
!                      vsstore(i,j,k,ispec)
           
!                 ! Stacey, a completer par la suite  
!                 rho_vp(i,j,k,ispec) = rhostore(i,j,k,ispec)*vpstore(i,j,k,ispec)
!                 rho_vs(i,j,k,ispec) = rhostore(i,j,k,ispec)*vsstore(i,j,k,ispec)
!                 !end pll
                
!                 !      kappastore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(ispec))* &
!                 !       (materials_ext_mesh(2,mat_ext_mesh(ispec))*materials_ext_mesh(2,mat_ext_mesh(ispec)) - &
!                 !        4.d0*materials_ext_mesh(3,mat_ext_mesh(ispec))*materials_ext_mesh(3,mat_ext_mesh(ispec))/3.d0)
!                 !      mustore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(ispec))*materials_ext_mesh(3,mat_ext_mesh(ispec))*&
!                 !  x    materials_ext_mesh(3,mat_ext_mesh(ispec))
!              enddo
!           enddo
!        enddo
!     enddo

  locval = 0
  ifseg = .false.
  xp = 0.d0
  yp = 0.d0
  zp = 0.d0

  do ispec=1,nspec
  ieoff = NGLLX * NGLLY * NGLLZ * (ispec-1)
  ilocnum = 0
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        ilocnum = ilocnum + 1
        xp(ilocnum+ieoff) = xstore(i,j,k,ispec)
        yp(ilocnum+ieoff) = ystore(i,j,k,ispec)
        zp(ilocnum+ieoff) = zstore(i,j,k,ispec)
      enddo
    enddo
  enddo
  enddo

  call get_global(nspec,xp,yp,zp,ibool,locval,ifseg,nglob,npointot, &
       minval(nodes_coords_ext_mesh(1,:)),maxval(nodes_coords_ext_mesh(1,:)))

  deallocate(xp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(yp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(zp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(locval,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(ifseg,stat=ier); if(ier /= 0) stop 'error in deallocate'

!
!- we can create a new indirect addressing to reduce cache misses
!
  allocate(copy_ibool_ori(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if(ier /= 0) stop 'error in allocate'
  allocate(mask_ibool(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'

  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:,:) = ibool(:,:,:,:)

  inumber = 0
  do ispec=1,nspec
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        if(mask_ibool(copy_ibool_ori(i,j,k,ispec)) == -1) then
! create a new point
          inumber = inumber + 1
          ibool(i,j,k,ispec) = inumber
          mask_ibool(copy_ibool_ori(i,j,k,ispec)) = inumber
        else
! use an existing point created previously
          ibool(i,j,k,ispec) = mask_ibool(copy_ibool_ori(i,j,k,ispec))
        endif
      enddo
    enddo
  enddo
  enddo

  deallocate(copy_ibool_ori,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(mask_ibool,stat=ier); if(ier /= 0) stop 'error in deallocate'

  allocate(xstore_dummy(nglob))
  allocate(ystore_dummy(nglob))
  allocate(zstore_dummy(nglob))
  do ispec = 1, nspec
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              iglobnum = ibool(i,j,k,ispec)
              xstore_dummy(iglobnum) = xstore(i,j,k,ispec)
           enddo
        enddo
     enddo
  enddo

  do ispec = 1, nspec
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              iglobnum = ibool(i,j,k,ispec)
              ystore_dummy(iglobnum) = ystore(i,j,k,ispec)
           enddo
        enddo
     enddo
  enddo

  do ispec = 1, nspec
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              iglobnum = ibool(i,j,k,ispec)
              zstore_dummy(iglobnum) = zstore(i,j,k,ispec)
           enddo
        enddo
     enddo
  enddo

! creating mass matrix (will be fully assembled with MPI in the solver)
  allocate(rmass(nglob))
  rmass(:) = 0._CUSTOM_REAL

  do ispec=1,nspec
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        weight=wxgll(i)*wygll(j)*wzgll(k)
        iglobnum=ibool(i,j,k,ispec)

        jacobianl=jacobianstore(i,j,k,ispec)

! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      rmass(iglobnum) = rmass(iglobnum) + &
          sngl((dble(rhostore(i,j,k,ispec)))  * dble(jacobianl) * weight)
    else
       rmass(iglobnum) = rmass(iglobnum) + rhostore(i,j,k,ispec) * jacobianl * weight
    endif

      enddo
    enddo
  enddo
  enddo  

 
  iboun(:,:) = .false. 
  do ispec2D = 1, nspec2D_xmin 
     iboun(1,ibelm_xmin(ispec2D)) = .true. 
  end do
  do ispec2D = 1, nspec2D_xmax 
     iboun(2,ibelm_xmax(ispec2D)) = .true. 
  end do
  do ispec2D = 1, nspec2D_ymin 
     iboun(3,ibelm_ymin(ispec2D)) = .true. 
  end do
  do ispec2D = 1, nspec2D_ymax 
     iboun(4,ibelm_ymax(ispec2D)) = .true. 
  end do
  do ispec2D = 1, NSPEC2D_BOTTOM
     iboun(5,ibelm_bottom(ispec2D)) = .true. 
  end do
  do ispec2D = 1, NSPEC2D_TOP
     iboun(6,ibelm_top(ispec2D)) = .true. 
  end do

  call get_jacobian_boundaries(myrank,iboun,nspec,xstore,ystore,zstore, &
       dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
       ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
       nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
       jacobian2D_xmin,jacobian2D_xmax, &
       jacobian2D_ymin,jacobian2D_ymax, &
       jacobian2D_bottom,jacobian2D_top, &
       normal_xmin,normal_xmax, &
       normal_ymin,normal_ymax, &
       normal_bottom,normal_top, &
       NSPEC2D_BOTTOM,NSPEC2D_TOP, &
       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

  call prepare_assemble_MPI (nelmnts_ext_mesh,ibool, &
       elmnts_ext_mesh, ESIZE, &
       nglob, &
       ninterface_ext_mesh, max_interface_size_ext_mesh, &
       my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
       ibool_interfaces_ext_mesh, &
       nibool_interfaces_ext_mesh &
       )

  ! Stacey put back
  call get_absorb_ext_mesh(myrank,iboun,nspec, &
       nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM)

! sort ibool comm buffers lexicographically
  allocate(nibool_interfaces_ext_mesh_true(ninterface_ext_mesh))

  do iinterface = 1, ninterface_ext_mesh

    allocate(xp(nibool_interfaces_ext_mesh(iinterface)))
    allocate(yp(nibool_interfaces_ext_mesh(iinterface)))
    allocate(zp(nibool_interfaces_ext_mesh(iinterface)))
    allocate(locval(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ifseg(nibool_interfaces_ext_mesh(iinterface)))
    allocate(reorder_interface_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ibool_interface_ext_mesh_dummy(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ind_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ninseg_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(iwork_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(work_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))

    do ilocnum = 1, nibool_interfaces_ext_mesh(iinterface)
      xp(ilocnum) = xstore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
      yp(ilocnum) = ystore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
      zp(ilocnum) = zstore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
    enddo

    call sort_array_coordinates(nibool_interfaces_ext_mesh(iinterface),xp,yp,zp, &
         ibool_interfaces_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
         reorder_interface_ext_mesh,locval,ifseg,nibool_interfaces_ext_mesh_true(iinterface), &
         ind_ext_mesh,ninseg_ext_mesh,iwork_ext_mesh,work_ext_mesh)

    deallocate(xp)
    deallocate(yp)
    deallocate(zp)
    deallocate(locval)
    deallocate(ifseg)
    deallocate(reorder_interface_ext_mesh)
    deallocate(ibool_interface_ext_mesh_dummy)
    deallocate(ind_ext_mesh)
    deallocate(ninseg_ext_mesh)
    deallocate(iwork_ext_mesh)
    deallocate(work_ext_mesh)

  enddo

! save the binary files
  call create_name_database(prname,myrank,LOCAL_PATH)
  open(unit=IOUT,file=prname(1:len_trim(prname))//'external_mesh.bin',status='unknown',action='write',form='unformatted')
  write(IOUT) nspec
  write(IOUT) nglob

  write(IOUT) xixstore
  write(IOUT) xiystore
  write(IOUT) xizstore
  write(IOUT) etaxstore
  write(IOUT) etaystore
  write(IOUT) etazstore
  write(IOUT) gammaxstore
  write(IOUT) gammaystore
  write(IOUT) gammazstore

  write(IOUT) jacobianstore

  !pll Stacey 
  write(IOUT) rho_vp
  write(IOUT) rho_vs
  write(IOUT) iflag_attenuation_store
  write(IOUT) NSPEC2DMAX_XMIN_XMAX 
  write(IOUT) NSPEC2DMAX_YMIN_YMAX
  write(IOUT) nimin
  write(IOUT) nimax
  write(IOUT) njmin
  write(IOUT) njmax
  write(IOUT) nkmin_xi 
  write(IOUT) nkmin_eta
  !end pll

  write(IOUT) kappastore
  write(IOUT) mustore

  write(IOUT) rmass

  write(IOUT) ibool

  write(IOUT) xstore_dummy
  write(IOUT) ystore_dummy
  write(IOUT) zstore_dummy

! boundary parameters
  write(IOUT) nspec2D_xmin
  write(IOUT) nspec2D_xmax
  write(IOUT) nspec2D_ymin
  write(IOUT) nspec2D_ymax
  write(IOUT) NSPEC2D_BOTTOM
  write(IOUT) NSPEC2D_TOP

  write(IOUT) ibelm_xmin
  write(IOUT) ibelm_xmax
  write(IOUT) ibelm_ymin
  write(IOUT) ibelm_ymax
  write(IOUT) ibelm_bottom
  write(IOUT) ibelm_top

  write(IOUT) normal_xmin
  write(IOUT) normal_xmax
  write(IOUT) normal_ymin
  write(IOUT) normal_ymax
  write(IOUT) normal_bottom
  write(IOUT) normal_top

  write(IOUT) jacobian2D_xmin
  write(IOUT) jacobian2D_xmax
  write(IOUT) jacobian2D_ymin
  write(IOUT) jacobian2D_ymax
  write(IOUT) jacobian2D_bottom
  write(IOUT) jacobian2D_top

! end boundary parameters

  write(IOUT) ninterface_ext_mesh
  write(IOUT) maxval(nibool_interfaces_ext_mesh)
  write(IOUT) my_neighbours_ext_mesh
  write(IOUT) nibool_interfaces_ext_mesh
  allocate(ibool_interfaces_ext_mesh_dummy(maxval(nibool_interfaces_ext_mesh),ninterface_ext_mesh))
  do i = 1, ninterface_ext_mesh
     ibool_interfaces_ext_mesh_dummy = ibool_interfaces_ext_mesh(1:maxval(nibool_interfaces_ext_mesh),:)
  enddo
  write(IOUT) ibool_interfaces_ext_mesh_dummy
  close(IOUT)

! compute the approximate amount of static memory needed to run the solver
  call memory_eval(nspec,nglob,maxval(nibool_interfaces_ext_mesh),ninterface_ext_mesh,static_memory_size)
  call max_all_dp(static_memory_size, max_static_memory_size)

  end subroutine create_regions_mesh_ext_mesh


!
!----
!

subroutine prepare_assemble_MPI (nelmnts,ibool, &
     knods, ngnode, &
     npoin, &
     ninterface, max_interface_size, &
     my_nelmnts_neighbours, my_interfaces, &
     ibool_interfaces_asteroid, &
     nibool_interfaces_asteroid &
     )

  implicit none

  include 'constants.h'

  integer, intent(in)  :: nelmnts, npoin, ngnode
  integer, dimension(ngnode,nelmnts), intent(in)  :: knods
  integer, dimension(NGLLX,NGLLY,NGLLZ,nelmnts), intent(in)  :: ibool

  integer  :: ninterface
  integer  :: max_interface_size
  integer, dimension(ninterface)  :: my_nelmnts_neighbours
  integer, dimension(6,max_interface_size,ninterface)  :: my_interfaces
  integer, dimension(NGLLX*NGLLX*max_interface_size,ninterface)  :: &
       ibool_interfaces_asteroid
  integer, dimension(ninterface)  :: &
       nibool_interfaces_asteroid

  integer  :: num_interface
  integer  :: ispec_interface

  logical, dimension(npoin)  :: mask_ibool_asteroid

  integer  :: ixmin, ixmax
  integer  :: iymin, iymax
  integer  :: izmin, izmax
  integer, dimension(ngnode)  :: n
  integer  :: e1, e2, e3, e4
  integer  :: type
  integer  :: ispec

  integer  :: k
  integer  :: npoin_interface_asteroid

  integer  :: ix,iy,iz


  ibool_interfaces_asteroid(:,:) = 0
  nibool_interfaces_asteroid(:) = 0

  do num_interface = 1, ninterface
     npoin_interface_asteroid = 0
     mask_ibool_asteroid(:) = .false.

     do ispec_interface = 1, my_nelmnts_neighbours(num_interface)
        ispec = my_interfaces(1,ispec_interface,num_interface)
        type = my_interfaces(2,ispec_interface,num_interface)
        do k = 1, ngnode
           n(k) = knods(k,ispec)
        end do
        e1 = my_interfaces(3,ispec_interface,num_interface)
        e2 = my_interfaces(4,ispec_interface,num_interface)
        e3 = my_interfaces(5,ispec_interface,num_interface)
        e4 = my_interfaces(6,ispec_interface,num_interface)
        call get_edge(ngnode, n, type, e1, e2, e3, e4, ixmin, ixmax, iymin, iymax, izmin, izmax)

        do iz = min(izmin,izmax), max(izmin,izmax)
           do iy = min(iymin,iymax), max(iymin,iymax)
              do ix = min(ixmin,ixmax), max(ixmin,ixmax)

                 if(.not. mask_ibool_asteroid(ibool(ix,iy,iz,ispec))) then
                    mask_ibool_asteroid(ibool(ix,iy,iz,ispec)) = .true.
                    npoin_interface_asteroid = npoin_interface_asteroid + 1
                    ibool_interfaces_asteroid(npoin_interface_asteroid,num_interface)=&
                         ibool(ix,iy,iz,ispec)
                 end if
              end do
           end do
        end do

     end do
     nibool_interfaces_asteroid(num_interface) = npoin_interface_asteroid


  end do

end subroutine prepare_assemble_MPI

!
!----
!

subroutine get_edge ( ngnode, n, type, e1, e2, e3, e4, ixmin, ixmax, iymin, iymax, izmin, izmax )

  implicit none

  include "constants.h"

  integer, intent(in)  :: ngnode
  integer, dimension(ngnode), intent(in)  :: n
  integer, intent(in)  :: type, e1, e2, e3, e4
  integer, intent(out)  :: ixmin, ixmax, iymin, iymax, izmin, izmax

  integer, dimension(4) :: en
  integer :: valence, i

   if ( type == 1 ) then
     if ( e1 == n(1) ) then
        ixmin = 1
        ixmax = 1
        iymin = 1
        iymax = 1
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(2) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = 1
        iymax = 1
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(3) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = NGLLY
        iymax = NGLLY
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(4) ) then
        ixmin = 1
        ixmax = 1
        iymin = NGLLY
        iymax = NGLLY
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(5) ) then
        ixmin = 1
        ixmax = 1
        iymin = 1
        iymax = 1
        izmin = NGLLZ
        izmax = NGLLZ
     end if
     if ( e1 == n(6) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = 1
        iymax = 1
        izmin = NGLLZ
        izmax = NGLLZ
     end if
     if ( e1 == n(7) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = NGLLY
        iymax = NGLLY
        izmin = NGLLZ
        izmax = NGLLZ
     end if
     if ( e1 == n(8) ) then
        ixmin = 1
        ixmax = 1
        iymin = NGLLY
        iymax = NGLLY
        izmin = NGLLZ
        izmax = NGLLZ
     end if
  else
     if ( type == 2 ) then
        if ( e1 ==  n(1) ) then
           ixmin = 1
           iymin = 1
           izmin = 1
           if ( e2 == n(2) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(4) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(5) ) then
              ixmax = 1
              iymax = 1
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(2) ) then
           ixmin = NGLLX
           iymin = 1
           izmin = 1
           if ( e2 == n(3) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(1) ) then
              ixmax = 1
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(6) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           end if

        end if
        if ( e1 == n(3) ) then
           ixmin = NGLLX
           iymin = NGLLY
           izmin = 1
           if ( e2 == n(4) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(2) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(7) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(4) ) then
           ixmin = 1
           iymin = NGLLY
           izmin = 1
           if ( e2 == n(1) ) then
              ixmax = 1
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(3) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(8) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(5) ) then
           ixmin = 1
           iymin = 1
           izmin = NGLLZ
           if ( e2 == n(1) ) then
              ixmax = 1
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(6) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           end if
           if ( e2 == n(8) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(6) ) then
           ixmin = NGLLX
           iymin = 1
           izmin = NGLLZ
           if ( e2 == n(2) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(7) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           end if
           if ( e2 == n(5) ) then
              ixmax = 1
              iymax = 1
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(7) ) then
           ixmin = NGLLX
           iymin = NGLLY
           izmin = NGLLZ
           if ( e2 == n(3) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(8) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           end if
           if ( e2 == n(6) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(8) ) then
           ixmin = 1
           iymin = NGLLY
           izmin = NGLLZ
           if ( e2 == n(4) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(5) ) then
              ixmax = 1
              iymax = 1
              izmax = NGLLZ
           end if
           if ( e2 == n(7) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if

     else
        if (type == 4) then
           en(1) = e1
           en(2) = e2
           en(3) = e3
           en(4) = e4

           valence = 0
           do i = 1, 4
              if ( en(i) == n(1)) then
                 valence = valence+1
              endif
              if ( en(i) == n(2)) then
                 valence = valence+1
              endif
              if ( en(i) == n(3)) then
                 valence = valence+1
              endif
              if ( en(i) == n(4)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = 1
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(1)) then
                 valence = valence+1
              endif
              if ( en(i) == n(2)) then
                 valence = valence+1
              endif
              if ( en(i) == n(5)) then
                 valence = valence+1
              endif
              if ( en(i) == n(6)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = 1
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(2)) then
                 valence = valence+1
              endif
              if ( en(i) == n(3)) then
                 valence = valence+1
              endif
              if ( en(i) == n(6)) then
                 valence = valence+1
              endif
              if ( en(i) == n(7)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = NGLLX
              iymin = 1
              izmin = 1
              ixmax = NGLLX
              iymax = NGLLZ
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(3)) then
                 valence = valence+1
              endif
              if ( en(i) == n(4)) then
                 valence = valence+1
              endif
              if ( en(i) == n(7)) then
                 valence = valence+1
              endif
              if ( en(i) == n(8)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = NGLLY
              izmin = 1
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(1)) then
                 valence = valence+1
              endif
              if ( en(i) == n(4)) then
                 valence = valence+1
              endif
              if ( en(i) == n(5)) then
                 valence = valence+1
              endif
              if ( en(i) == n(8)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = 1
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(5)) then
                 valence = valence+1
              endif
              if ( en(i) == n(6)) then
                 valence = valence+1
              endif
              if ( en(i) == n(7)) then
                 valence = valence+1
              endif
              if ( en(i) == n(8)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = NGLLZ
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           endif

        else
           stop 'ERROR get_edge'
        endif

     end if
  end if

end subroutine get_edge



!pll
! subroutine interface(iflag,flag_below,flag_above,ispec,nspec,i,j,k,xstore,ystore,zstore,ibedrock)

! implicit none

! include "constants.h"

! integer :: iflag,flag_below,flag_above
! integer :: ispec,nspec
! integer :: i,j,k
! double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore
! real(kind=CUSTOM_REAL), dimension(NX_TOPO_ANT,NY_TOPO_ANT) :: ibedrock
! integer, parameter :: NUMBER_OF_STATIONS = 1
! double precision, parameter :: RADIUS_TO_EXCLUDE = 250.d0
! double precision, dimension(NUMBER_OF_STATIONS) :: utm_x_station,utm_y_station

! !-------------------

! !for Piero
! logical :: is_around_a_station
! integer :: istation

! ! store bedrock values
! integer ::  icornerlat,icornerlong
! double precision ::  lat,long,elevation_bedrock
! double precision ::  lat_corner,long_corner,ratio_xi,ratio_eta


! !! DK DK store the position of the six stations to be able to
! !! DK DK exclude circles around each station to make sure they are on the bedrock
! !! DK DK and not in the ice
!    utm_x_station(1) =  783500.6250000d0
!    utm_y_station(1) = -11828.7519531d0

!    utm_x_station(2) =  853644.5000000d0
!    utm_y_station(2) = -114.0138092d0

!    utm_x_station(3) = 863406.0000000d0
!    utm_y_station(3) = -53736.1640625d0

!    utm_x_station(4) =   823398.8125000d0
!    utm_y_station(4) = 29847.4511719d0

!    utm_x_station(5) = 863545.3750000d0
!    utm_y_station(5) = 19669.6621094d0

!    utm_x_station(6) =  817099.3750000d0
!    utm_y_station(6) = -24430.2871094d0

! ! since we have suppressed UTM projection for Piero Basini, UTMx is the same as long
! ! and UTMy is the same as lat
!     long = xstore(i,j,k,ispec)
!     lat =  ystore(i,j,k,ispec)

! ! get coordinate of corner in model
!     icornerlong = int((long - ORIG_LONG_TOPO_ANT) / DEGREES_PER_CELL_TOPO_ANT) + 1
!     icornerlat = int((lat - ORIG_LAT_TOPO_ANT) / DEGREES_PER_CELL_TOPO_ANT) + 1

! ! avoid edge effects and extend with identical point if outside model
!     if(icornerlong < 1) icornerlong = 1
!     if(icornerlong > NX_TOPO_ANT-1) icornerlong = NX_TOPO_ANT-1
!     if(icornerlat < 1) icornerlat = 1
!     if(icornerlat > NY_TOPO_ANT-1) icornerlat = NY_TOPO_ANT-1

! ! compute coordinates of corner
!     long_corner = ORIG_LONG_TOPO_ANT + (icornerlong-1)*DEGREES_PER_CELL_TOPO_ANT
!     lat_corner = ORIG_LAT_TOPO_ANT + (icornerlat-1)*DEGREES_PER_CELL_TOPO_ANT

! ! compute ratio for interpolation
!     ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO_ANT
!     ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO_ANT

! ! avoid edge effects
!     if(ratio_xi < 0.) ratio_xi = 0.
!     if(ratio_xi > 1.) ratio_xi = 1.
!     if(ratio_eta < 0.) ratio_eta = 0.
!     if(ratio_eta > 1.) ratio_eta = 1.

! ! interpolate elevation at current point
!     elevation_bedrock = &
!       ibedrock(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
!       ibedrock(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
!       ibedrock(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
!       ibedrock(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

! !! DK DK exclude circles around each station to make sure they are on the bedrock
! !! DK DK and not in the ice
!   is_around_a_station = .false.
!   do istation = 1,NUMBER_OF_STATIONS
!     if(sqrt((xstore(i,j,k,ispec) - utm_x_station(istation))**2 + (ystore(i,j,k,ispec) - &
!          utm_y_station(istation))**2) < RADIUS_TO_EXCLUDE) then
!       is_around_a_station = .true.
!       exit
!     endif
!   enddo

! ! we are above the bedrock interface i.e. in the ice, and not too close to a station
!   if(zstore(i,j,k,ispec) >= elevation_bedrock .and. .not. is_around_a_station) then
!      iflag = flag_above
!      !iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_ICE
!      ! we are below the bedrock interface i.e. in the bedrock, or close to a station
!   else
!      iflag = flag_below
!      !iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_BEDROCK
!   endif
    

! end subroutine interface
