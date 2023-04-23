!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage < http://www.axisem.info>
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
!    along with AxiSEM.  If not, see < http://www.gnu.org/licenses/>.
!

!=========================================================================================
module data_mesh

  ! Arrays here pertain to some sort of mesh peculiarities and mainly serve
  ! as information or parameters for many "if"-decisions such as
  ! axis, north, element type, solid-fluid boundary mapping, coarsening,
  ! and specifically also related to the background model such as
  ! solid-fluid boundary mapping, discontinuities, and element arrays to
  ! map solid/fluid to global element domains.
  ! These quantities are in active memory throughout the simulation.
  !
  ! Any global arrays containing properties inside elements are defined in data_matr.

  use global_parameters
  implicit none

  public

  ! Very basic mesh parameters, have been in mesh_params.h before
  integer , protected ::         npol ! < polynomial order
  integer , protected ::        nelem ! < proc. els
  integer , protected ::       npoint ! < proc. all pts
  integer , protected ::    nel_solid ! < proc. solid els
  integer , protected ::    nel_fluid ! < proc. fluid els
  integer , protected :: npoint_solid ! < proc. solid pts
  integer , protected :: npoint_fluid ! < proc. fluid pts
  integer , protected ::  nglob_solid ! < proc. slocal pts
  integer , protected ::  nglob_fluid ! < proc. flocal pts
  integer , protected ::     nel_bdry ! < proc. solid-fluid bndry els
  integer , protected ::        ndisc ! < # disconts in bkgrd model
  integer , protected ::   nproc_mesh ! < number of processors
  integer , protected :: lfbkgrdmodel ! < length of bkgrdmodel name

  ! global number in solid varies across procs due to central cube domain decomposition
  integer                                       :: nglob
  ! global numbering array for the solid and fluid assembly
  integer, protected, allocatable, dimension(:) :: igloc_solid ! (npoint_solid)
  integer, protected, allocatable, dimension(:) :: igloc_fluid ! (npoint_fluid)

  ! Misc definitions
  integer                           :: nsize
  logical                           :: do_mesh_tests

  real(kind=realkind), allocatable  :: gvec_solid(:,:)
  real(kind=realkind), allocatable  :: gvec_fluid(:)

  ! Deprecated elemental mesh (radius & colatitude of elemental midpoint)
  ! This is used for blow up localization. Might want to remove this when
  ! everything is running smoothly in all eternities...
  real(kind=realkind), allocatable  :: mean_rad_colat_solid(:,:)
  real(kind=realkind), allocatable  :: mean_rad_colat_fluid(:,:)

  ! Global mesh informations
  real(kind=dp)                     :: router ! Outer radius (surface)

  ! critical mesh parameters (spacing/velocity, characteristic lead time etc)
  real(kind=dp)                     :: pts_wavelngth
  real(kind=dp)                     :: hmin_glob, hmax_glob
  real(kind=dp)                     :: min_distance_dim, min_distance_nondim
  real(kind=dp)                     :: char_time_max
  integer                           :: char_time_max_globel
  real(kind=dp)                     :: char_time_max_rad, char_time_max_theta
  real(kind=dp)                     :: char_time_min
  integer                           :: char_time_min_globel
  real(kind=dp)                     :: char_time_min_rad, char_time_min_theta
  real(kind=dp)                     :: vpmin, vsmin, vpmax, vsmax
  real(kind=dp)                     :: vpminr, vsminr, vpmaxr, vsmaxr
  integer, dimension(3)             :: vpminloc, vsminloc, vpmaxloc, vsmaxloc

  !----------------------------------------------------------------------
  ! Axial elements
  integer                           :: naxel, naxel_solid, naxel_fluid
  integer, protected, allocatable   :: ax_el(:), ax_el_solid(:), ax_el_fluid(:)
  logical,            allocatable   :: axis_solid(:)
  logical,            allocatable   :: axis_fluid(:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Background Model related
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solid-Fluid boundary----------------------------------------------------

  ! mapping nel_bdry elements into the solid/fluid element numbers
  integer, protected, allocatable :: bdry_solid_el(:), bdry_fluid_el(:)

  ! mapping the z coordinate of the boundary for each element
  ! (depending on north/south, above/below)
  integer, protected, allocatable :: bdry_jpol_solid(:), bdry_jpol_fluid(:)

  ! Boolean to determine whether proc has solid-fluid boundary elements
  logical :: have_bdry_elem

  !Not used anywhere
  !! integer array of size nel_bdry containing the "global" element number
  !! for 1:nel_bdry
  !integer, dimension(nel_bdry) :: ibdryel


  ! Background model--------------------------------------------------------
  character(len=100)          :: bkgrdmodel
  character(len=100)          :: meshname
  logical                     :: have_fluid
  real(kind=dp), allocatable  :: discont(:)
  logical, allocatable        :: solid_domain(:)
  integer, allocatable        :: idom_fluid(:)
  real(kind=dp)               :: rmin, minh_ic, maxh_ic, maxh_icb
  logical                     :: anel_true ! anelastic model?
  !--------------------------------------------------------------------------

  ! Receiver locations
  integer                      :: maxind_glob, maxind, ind_first, ind_last
  integer                      :: num_rec, num_rec_tot
  integer, allocatable         :: surfelem(:), jsurfel(:)
  real(kind=sp), allocatable   :: surfcoord(:)
  integer                      :: ielepi, ielantipode, ielequ
  integer, allocatable         :: recfile_el(:,:), loc2globrec(:)
  logical                      :: have_epi, have_equ, have_antipode
  real                         :: dtheta_rec

  ! CMB receivers (same as receivers, just above CMB instead)
  integer                      :: num_cmb
  integer, allocatable         :: cmbfile_el(:,:), loc2globcmb(:)
  !--------------------------------------------------------------------------

  ! for xdmf plotting
  integer                      :: nelem_plot, npoint_plot
  logical, allocatable         :: plotting_mask(:,:,:)
  integer, allocatable         :: mapping_ijel_iplot(:,:,:)

  ! for kernel wavefields in displ_only mode
  integer                      :: nelem_kwf_global, nelem_kwf
  integer                      :: npoint_kwf_global, npoint_kwf
  integer                      :: npoint_solid_kwf, npoint_fluid_kwf
  logical, allocatable         :: kwf_mask(:,:,:)
  integer, allocatable         :: mapping_ijel_ikwf(:,:,:)
  integer, allocatable         :: midpoint_mesh_kwf(:), eltype_kwf(:), axis_kwf(:)
  integer, allocatable         :: fem_mesh_kwf(:,:)
  integer, allocatable         :: sem_mesh_kwf(:,:,:)

  ! Only needed before the simulation and later deallocated
  ! Global mesh informations
  integer, dimension(:,:), allocatable          :: lnods ! (nelem,1:8)
  character(len=6), dimension(:), allocatable   :: eltype ! (nelem)
  integer                                       :: npoin
  real(kind=dp), dimension(:,:), allocatable    :: crd_nodes ! (npoin,2)
  logical, dimension(:), allocatable            :: coarsing,north ! (nelem)
  integer                                       :: num_spher_radii
  real(kind=dp)   , dimension(:), allocatable   :: spher_radii
  ! Axial elements
  logical, dimension(:), allocatable            :: axis ! (nelem)

  ! Mapping between solid/fluid elements:
  ! integer array of size nel_fluid containing the glocal (global per-proc)
  ! element number for 1:nel_solid/fluid
  integer, dimension(:), allocatable            :: ielsolid ! (nel_solid)
  integer, dimension(:), allocatable            :: ielfluid ! (nel_fluid)

contains

!-----------------------------------------------------------------------------------------
!> Read parameters formerly in mesh_params.h
!! It is slightly dirty to have this routine in a data module
!! but it allows to define the variables as 'protected', i.e.
!! fixed outside of this module.
subroutine read_mesh_basics(iounit)

   integer, intent(in)   :: iounit

   read(iounit) nproc_mesh
   read(iounit) npol
   read(iounit) nelem
   read(iounit) npoint
   read(iounit) nel_solid
   read(iounit) nel_fluid
   read(iounit) npoint_solid
   read(iounit) npoint_fluid
   read(iounit) nglob_solid
   read(iounit) nglob_fluid
   read(iounit) nel_bdry
   read(iounit) ndisc
   read(iounit) lfbkgrdmodel

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_mesh_advanced(iounit)
   use data_io, only: verbose
   use data_spec
   integer, intent(in)  :: iounit
   integer              :: iptcp, iel, inode

   allocate(eta(0:npol))
   allocate(dxi(0:npol))
   allocate(wt(0:npol))
   allocate(xi_k(0:npol))
   allocate(wt_axial_k(0:npol))
   allocate(G1(0:npol,0:npol))
   allocate(G1T(0:npol,0:npol))
   allocate(G2(0:npol,0:npol))
   allocate(G2T(0:npol,0:npol))
   allocate(G0(0:npol))

   ! spectral stuff
   read(iounit) xi_k
   read(iounit) eta
   read(iounit) dxi
   read(iounit) wt
   read(iounit) wt_axial_k
   read(iounit) G0
   read(iounit) G1
   read(iounit) G1T
   read(iounit) G2
   read(iounit) G2T

   read(iounit) npoin

   if (verbose > 1) then
      write(69,*) 'reading coordinates/control points...'
      write(69,*) 'global number of control points:',npoin
   endif
   allocate(crd_nodes(1:npoin,1:2))

   read(iounit) crd_nodes(:,1)
   read(iounit) crd_nodes(:,2)
   do iptcp = 1, npoin
      if (abs(crd_nodes(iptcp,2)) < 1.e-8) crd_nodes(iptcp,2) = zero
   enddo

   allocate(lnods(1:nelem,1:8))
   do iel = 1, nelem
      read(iounit) (lnods(iel,inode), inode=1,8)
   enddo


   ! Number of global distinct points (slightly differs for each processor!)
   read(iounit) nglob
   if (verbose > 1) write(69,*) '  global number:', nglob

   ! Element type
   allocate(eltype(1:nelem), coarsing(1:nelem))
   read(iounit) eltype
   read(iounit) coarsing

   !!!!!!!!!!! SOLID/FLUID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! mapping from sol/flu (1:nel_fluid) to global element numbers (1:neltot)
   if (verbose > 1) write(69,*) 'reading solid/fluid domain info...'
   allocate(ielsolid(1:nel_solid))
   allocate(ielfluid(1:nel_fluid))
   read(iounit) ielsolid
   read(iounit) ielfluid

   ! slocal numbering
   allocate(igloc_solid(npoint_solid))
   read(iounit) igloc_solid(1:npoint_solid)

   ! flocal numbering
   allocate(igloc_fluid(npoint_fluid))
   read(iounit) igloc_fluid(1:npoint_fluid)

   ! Solid-Fluid boundary
   if (verbose > 1) write(69,*) 'reading solid/fluid boundary info...'
   read(iounit) have_bdry_elem

   if (have_bdry_elem) then
      allocate(bdry_solid_el(1:nel_bdry))
      allocate(bdry_fluid_el(1:nel_bdry))
      allocate(bdry_jpol_solid(1:nel_bdry))
      allocate(bdry_jpol_fluid(1:nel_bdry))
      read(iounit) bdry_solid_el(1:nel_bdry)
      read(iounit) bdry_fluid_el(1:nel_bdry)
      read(iounit) bdry_jpol_solid(1:nel_bdry)
      read(iounit) bdry_jpol_fluid(1:nel_bdry)
   endif
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_mesh_axel(iounit)
   integer, intent(in) :: iounit

   allocate(ax_el(naxel))
   allocate(ax_el_solid(1:naxel_solid))
   allocate(ax_el_fluid(1:naxel_fluid))
   allocate(axis_solid(nel_solid))
   allocate(axis_fluid(nel_fluid))

   read(iounit) ax_el(1:naxel)
   read(iounit) ax_el_solid(1:naxel_solid)
   read(iounit) ax_el_fluid(1:naxel_fluid)
end subroutine
!-----------------------------------------------------------------------------------------

end module data_mesh
!=========================================================================================
