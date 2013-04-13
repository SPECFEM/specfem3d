!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
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

! interface model file
! example file only, unused so far

! !  Piero
!  module bedrock
!
!  real,dimension(:,:),allocatable :: ibedrock
!
!  end module bedrock

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine model_bedrock_broadcast(myrank)
!
!! standard routine to setup model
!
!  use bedrock
!
!  implicit none
!
!  include "constants.h"
!  ! standard include of the MPI library
!  include 'mpif.h'
!
!  integer :: myrank
!
!  ! local parameters
!  integer :: idummy
!
!  ! dummy to ignore compiler warnings
!  idummy = myrank
!
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

!  end subroutine model_bedrock_broadcast

!
!-------------------------------------------------------------------------------------------------
!

!
!  subroutine read_bedrock_file()
!
!  use bedrock
!
!  implicit none
!
!  include "constants.h"
!---
!
! ADD YOUR MODEL HERE
!
!---
!
!  end subroutine read_external_model


!
!-------------------------------------------------------------------------------------------------
!


!  subroutine model_bedrock_store()
!
! use bedrock
!
! implicit none
!
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


!            if(idoubling(ispec) == IFLAG_ONE_LAYER_TOPOGRAPHY .or. &
!               idoubling(ispec) == IFLAG_BEDROCK_down_to_14km) then

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
!                  qmu_attenuation_store(i,j,k,ispec) = 1.0 ! IATTENUATION_ICE

!                  ! we are below the bedrock interface i.e. in the bedrock, or close to a station
!               else
!                  vp = 5800.d0
!                  vs = 3200.d0
!                  rho = 2600.d0
!                  qmu_attenuation_store(i,j,k,ispec) = 9000.0 ! IATTENUATION_BEDROCK
!               endif

!            else if(idoubling(ispec) == IFLAG_14km_to_23p4km) then
!               vp = 6800.d0
!               vs = 3900.d0
!               rho = 2900.d0
!               qmu_attenuation_store(i,j,k,ispec) = 9000.0 ! IATTENUATION_BEDROCK

!            else if(idoubling(ispec) == IFLAG_MANTLE_BELOW_23p4km) then
!               vp = 8100.d0
!               vs = 4480.d0
!               rho = 3380.d0
!               qmu_attenuation_store(i,j,k,ispec) = 9000.0 ! IATTENUATION_BEDROCK

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
!                     endif

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
!                 !      mustore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(ispec))* &
!                                                         materials_ext_mesh(3,mat_ext_mesh(ispec))*&
!                 !  x    materials_ext_mesh(3,mat_ext_mesh(ispec))
!              enddo
!           enddo
!        enddo
!     enddo
!
!  end subroutine


!
!-------------------------------------------------------------------------------------------------
!


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
!      !qmu_attenuation_store(i,j,k,ispec) = 1.0 ! IATTENUATION_ICE
!      ! we are below the bedrock interface i.e. in the bedrock, or close to a station
!   else
!      iflag = flag_below
!      !qmu_attenuation_store(i,j,k,ispec) = 9000.0 ! IATTENUATION_BEDROCK
!   endif


! end subroutine interface
