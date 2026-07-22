!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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


module decompose_mesh_par

  use constants, only: IIN_DB,MAX_STRING_LEN,NGNOD_EIGHT_CORNERS

  use shared_parameters

  use part_decompose_mesh, only: ACOUSTIC_LOAD, &
    write_interfaces_database,write_moho_surface_database,write_glob2loc_nodes_database, &
    write_material_props_database,write_boundaries_database, &
    write_partition_database,write_cpml_database, &
    acoustic_elastic_poro_load,mesh2dual_ncommonnodes, &
    build_glob2loc_elmnts,build_glob2loc_nodes,build_interfaces,poro_elastic_repartitioning,moho_surface_repartitioning

  !! setting up wavefield discontinuity interface
  use part_decompose_mesh, only: write_wavefield_discontinuity_database
  implicit none

! note: the poroelastic repartitioning routine to parallelize the poroelastic-elastic interface
!       might break the load balancing created by the domain decomposer for high-performance computing
!
!       however, poroelastic-elastic interface need to have coupled poroelastic and elastic elements in the same
!       partition in order to compute proper coupling terms for now...
  logical, parameter :: PORO_INTERFACE_REPARTITIONING = .true.

  ! number of partitions
  integer :: nparts = 0

  ! mesh arrays
  integer :: nspec = 0
  integer, dimension(:,:), allocatable  :: elmnts
  integer, dimension(:,:), allocatable  :: mat
  integer, dimension(:), allocatable  :: part

  integer :: nnodes = 0
  double precision, dimension(:,:), allocatable  :: nodes_coords

  integer, dimension(:), allocatable  :: xadj
  integer, dimension(:), allocatable  :: adjncy
  integer, dimension(:), allocatable  :: nnodes_elmnts
  integer, dimension(:), allocatable  :: nodes_elmnts
  integer, dimension(:), allocatable  :: elmnts_load

  integer, dimension(:), pointer  :: glob2loc_elmnts
  integer, dimension(:), pointer  :: glob2loc_nodes_nparts
  integer, dimension(:), pointer  :: glob2loc_nodes_parts
  integer, dimension(:), pointer  :: glob2loc_nodes

  integer, dimension(:), pointer  :: tab_size_interfaces, tab_interfaces
  integer, dimension(:), allocatable  :: my_interfaces
  integer, dimension(:), allocatable  :: my_nb_interfaces
  integer :: ninterfaces = 0

  integer :: nsize           ! max number of elements that contain the same node
  integer :: nb_edges

  integer :: sup_neighbor   ! majoration (overestimate) of the maximum number of neighbors per element

  integer :: ipart, nnodes_loc, nspec_local, ncommonnodes
  integer :: num_elmnt, num_node, num_mat

  ! boundaries
  integer :: nspec2D_xmin = 0
  integer :: nspec2D_xmax = 0
  integer :: nspec2D_ymin = 0
  integer :: nspec2D_ymax = 0
  integer :: nspec2D_bottom = 0
  integer :: nspec2D_top = 0
  integer, dimension(:), allocatable :: ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top
  integer, dimension(:,:), allocatable :: nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin
  integer, dimension(:,:), allocatable :: nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top

  ! C-PML absorbing boundary conditions
  integer :: nspec_cpml = 0
  integer, dimension(:), allocatable :: CPML_to_spec, CPML_regions
  logical, dimension(:), allocatable :: is_CPML

  ! moho surface (optional)
  integer :: nspec2D_moho = 0
  integer, dimension(:), allocatable :: ibelm_moho
  integer, dimension(:,:), allocatable :: nodes_ibelm_moho

  character(len=MAX_STRING_LEN) :: prname

  !pll
  double precision , dimension(:,:), allocatable :: mat_prop
  integer :: count_def_mat,count_undef_mat,imat
  character(len=MAX_STRING_LEN), dimension(:,:), allocatable :: undef_mat_prop

  ! default mesh file directory
  character(len=MAX_STRING_LEN) :: localpath_name
  character(len=MAX_STRING_LEN) :: outputpath_name

  integer, parameter :: IIN_database = IIN_DB

  ! LTS simulations
  ! element p-refinement values (like 1 2 4 8 ..; p == 1 being coarsest, p == 8 finer local time step dt/p )
  integer,dimension(:),allocatable :: ispec_p_refine
  integer,dimension(:),allocatable :: p_level
  integer,dimension(:),allocatable :: num_ispec_level
  integer :: num_p_level = 0

  ! wavefield discontinuity
  !! boundary of wavefield discontinuity, read from database file
  integer :: nb_wd

  !! boundary_to_ispec_wd(nb_wd)
  !! the element the boundary belongs to, read from database file
  !! each point on the boundary belongs to two sides of the boundary
  !! here the element must be on the inner side of the boundary
  integer, dimension(:), allocatable :: boundary_to_ispec_wd

  !! side_wd(nb_wd)
  !! integers specifying which side the boundary is in the element
  !! read from database file
  !! side_wd = 1--8: only one vertex is on the boundary
  !! side_wd = 9--20: only one edge is on the boundary
  !! side_wd = 21--26: one face is on the boundary
  integer, dimension(:), allocatable :: side_wd

end module decompose_mesh_par

