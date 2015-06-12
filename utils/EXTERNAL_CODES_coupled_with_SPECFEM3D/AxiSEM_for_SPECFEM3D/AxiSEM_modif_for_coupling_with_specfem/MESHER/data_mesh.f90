!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

module data_mesh

  use global_parameters
  
  implicit none
  public

  !skeleton parameters
  integer :: neltot
  integer :: npointot
  real(kind=dp)   , dimension(:), allocatable   :: sg,zg
  integer, dimension(:,:), allocatable          :: lnodesg
  character(len=6), dimension(:), allocatable   :: eltypeg
  logical, dimension(:), allocatable            :: coarsing
  
  !Region to which the element pertains, in the case of a stratified bkgrd model
  integer, dimension(:), allocatable            :: region, nel_region
  real(kind=dp)   , dimension(:), allocatable   :: scom, zcom, thetacom, rcom 
                                                   ! com = center of mass
  
  ! axial elements
  integer, allocatable  :: naxelp(:), naxel_solidp(:), naxel_fluidp(:)
  integer, allocatable  :: ax_elp(:,:), ax_el_solidp(:,:), ax_el_fluidp(:,:)
  integer,allocatable   :: axis(:,:), axis_solid(:,:), axis_fluid(:,:)
  
  ! Solid fluid distinction
  integer                            :: neltot_fluid
  logical, dimension(:), allocatable :: fluid
  integer, dimension(:), allocatable :: ielem_fluid, inv_ielem_fluid

  integer                            :: neltot_solid
  logical, dimension(:), allocatable :: solid
  integer, dimension(:), allocatable :: ielem_solid, inv_ielem_solid
  
  !Boundary informations (ICB, CMB, from solid and fluid perspectives)
  integer                               :: nbcnd    ! Number of boundaries
  integer, dimension(:),allocatable     :: nbelem 
  integer, dimension(:,:), allocatable  :: belem    ! List of boundary elements
  integer, dimension(:,:),allocatable   :: my_neighbour
  
  ! Solid-fluid boundary arrays (needed by solver)
  integer, dimension(:,:),allocatable  :: bdry_above_el, bdry_below_el 
  integer, dimension(:,:),allocatable  :: bdry_solid_el, bdry_fluid_el 
  
  real(kind=dp)   , dimension(:,:,:), allocatable :: bdry_s,bdry_z
  integer, dimension(:,:,:),allocatable  :: bdry_globnum_above
  integer, dimension(:,:,:),allocatable  :: bdry_globnum_below
  integer, dimension(:,:,:),allocatable  :: bdry_locnum_above
  integer, dimension(:,:,:),allocatable  :: bdry_locnum_below
  
  integer, dimension(:,:,:),allocatable  :: bdry_globnum_solid
  integer, dimension(:,:,:),allocatable  :: bdry_globnum_fluid
  integer, dimension(:,:,:),allocatable  :: bdry_locnum_solid
  integer, dimension(:,:,:),allocatable  :: bdry_locnum_fluid
  
  integer, dimension(:,:), allocatable   :: bdry_jpol_solid, bdry_jpol_fluid
  
  ! central region
  real(kind=dp)   , allocatable :: s_arr(:,:), z_arr(:,:)
  integer, allocatable          :: central_is_iz_to_globiel(:,:)
contains

subroutine empty_data_mesh
  deallocate(rcom, scom, zcom, thetacom)
end subroutine empty_data_mesh

end module data_mesh
