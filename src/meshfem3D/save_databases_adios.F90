!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

!==============================================================================
!> \file save_databases_adios.F90
!!
!! \author MPBL
!==============================================================================
#include "config.fh"


!==============================================================================
subroutine save_databases_adios(LOCAL_PATH, myrank, sizeprocs, &
                                nspec,nglob,iproc_xi,iproc_eta, &
                                NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta, &
                                ibool,nodes_coords,ispec_material_id, &
                                nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                NSPEC2D_BOTTOM,NSPEC2D_TOP, NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                                ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                NMATERIALS,material_properties, &
                                nspec_CPML,CPML_to_spec,CPML_regions,is_CPML)

  use constants, only: MAX_STRING_LEN,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,ADIOS_TRANSPORT_METHOD, &
    NGLLX,NGLLY,NGLLZ

  use adios_helpers_mod, only: define_adios_global_array1d,define_adios_scalar, &
    write_adios_global_1d_array,write_adios_global_string_1d_array, &
    define_adios_local_string_1d_array

  use safe_alloc_mod, only: safe_alloc,safe_dealloc

  implicit none

  include "constants_meshfem3D.h"

  ! MPI variables
  integer :: myrank, sizeprocs

  ! number of spectral elements in each block
  integer nspec

  ! number of vertices in each block
  integer nglob

  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX),
  ! S for South (= ETA_MIN), N for North (= ETA_MAX)
  integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8
  integer iproc_xi,iproc_eta
  integer NPROC_XI,NPROC_ETA
  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)
  integer addressing(0:NPROC_XI-1,0:NPROC_ETA-1)

  ! arrays with the mesh
  integer ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)
  double precision :: nodes_coords(nglob,3)

  integer ispec_material_id(nspec)

  ! boundary parameters locator
  integer NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer ibelm_bottom(NSPEC2D_BOTTOM)
  integer ibelm_top(NSPEC2D_TOP)

  ! material properties
  integer :: NMATERIALS
  ! first dimension  : material_id
  ! second dimension : #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
  double precision , dimension(NMATERIALS,NUMBER_OF_MATERIAL_PROPERTIES) :: material_properties

  ! CPML
  integer, intent(in) :: nspec_CPML
  integer, dimension(nspec_CPML), intent(in) :: CPML_to_spec,CPML_regions
  logical, dimension(nspec), intent(in) :: is_CPML
  integer :: nspec_CPML_total

  integer :: i,ispec,ier

  ! name of the database files
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! for MPI interfaces
  integer ::  nb_interfaces,nspec_interfaces_max
  logical, dimension(8) ::  interfaces
  integer, dimension(8) ::  nspec_interface

  ! materials
  integer :: ndef,nundef
  integer :: icount,is,ie
  integer :: mat_id,domain_id
  integer(kind=4), dimension(2,nspec) :: material_index
  double precision , dimension(:,:),allocatable :: matpropl
  character(len=MAX_STRING_LEN*6*NMATERIALS) :: undef_matpropl

  integer :: ngnod, ngnod2d

  !--- Local parameters for ADIOS ---
  character(len=MAX_STRING_LEN) :: output_name
  character(len=*), parameter :: group_name = "SPECFEM3D_DATABASES"
  integer(kind=8) :: group, handle
  integer(kind=8) :: groupsize, totalsize
  integer :: local_dim

  !--- Variables to allreduce - wmax stands for world_max
  integer :: nglob_wmax, nspec_wmax, nmaterials_wmax, &
             nspec2d_xmin_wmax, nspec2d_xmax_wmax, &
             nspec2d_ymin_wmax, nspec2d_ymax_wmax, &
             nspec2d_bottom_wmax, nspec2d_top_wmax, &
             nb_interfaces_wmax, nspec_interfaces_max_wmax
  integer, parameter :: num_vars = 11
  integer, dimension(num_vars) :: max_global_values

  !--- Temporary arrays for writes
  integer, dimension(:,:), allocatable :: nodes_ibelm_xmin, nodes_ibelm_xmax, &
                                          nodes_ibelm_ymin, nodes_ibelm_ymax, &
                                          nodes_ibelm_bottom, nodes_ibelm_top
  integer, dimension(:), allocatable :: neighbors_mesh, num_elmnts_mesh
  integer, dimension(:,:,:), allocatable :: interfaces_mesh
  integer, dimension(:,:), allocatable :: elmnts_mesh
  integer :: interface_num, ispec_interface
  integer :: comm

  !---------------------------.
  ! Setup the values to write |
  !---------------------------'
  ngnod   = NGLLX_M * NGLLY_M *  NGLLZ_M
  ngnod2d = NGLLX_M * NGLLY_M

  ! assignes material index
  ! format: (1,ispec) = #material_id , (2,ispec) = #material_definition
  material_index (:,:) = 0
  do ispec = 1, nspec
    ! material id
    material_index(1,ispec) = ispec_material_id(ispec)
    ! material definition: 1 = interface type / 2 = tomography type
    if (ispec_material_id(ispec) > 0) then
      ! dummy value, not used any further
      material_index(2,ispec) = 1
    else
      ! negative material ids
      ! by default, assumes tomography model material
      ! (note: interface type not implemented yet...)
      material_index(2,ispec) = 2
    endif
  enddo

  ! Materials properties
  ! counts defined/undefined materials
  ndef = 0
  nundef = 0
  do i = 1,NMATERIALS
    mat_id = material_properties(i,8)
    if (mat_id > 0) ndef = ndef + 1
    if (mat_id < 0) nundef = nundef + 1
  enddo

  ! defined material properties
  ! pad dummy zeros to fill up 16 entries (poroelastic medium not allowed)
  call safe_alloc(matpropl,16, ndef,"matpropl")
  matpropl(:,:) = 0.d0
  icount = 0
  do i = 1,NMATERIALS
    mat_id = material_properties(i,8)
    if (mat_id > 0) then
      icount = icount + 1
      matpropl(1:7, icount) = material_properties(i,1:7)
    endif
  enddo
  if (icount /= ndef) stop 'Error icount not equal to ndef'

  ! undefined material properties
  ! we use a 1D character string for ADIOS read/write
  undef_matpropl(:) = ''
  icount = 0
  do i = 1,NMATERIALS
    domain_id = material_properties(i,7)
    mat_id = material_properties(i,8)
    if (mat_id < 0) then
      icount = icount + 1
      ! format:
      ! #material_id #type-keyword #domain-name #tomo-filename #tomo_id #domain_id
      ! format example tomography: -1 tomography elastic tomography_model.xyz 0 2
      ! material id
      is = (icount-1)*6*MAX_STRING_LEN + 1
      ie = is + MAX_STRING_LEN
      write(undef_matpropl(is:ie),*) mat_id
      ! name
      is = (icount-1)*6*MAX_STRING_LEN + (1*MAX_STRING_LEN + 1)
      ie = is + MAX_STRING_LEN
      undef_matpropl(is:ie) = 'tomography'
      ! name type
      is = (icount-1)*6*MAX_STRING_LEN + (2*MAX_STRING_LEN + 1)
      ie = is + MAX_STRING_LEN
      select case (domain_id)
      case (IDOMAIN_ACOUSTIC)
        undef_matpropl(is:ie) = 'acoustic'
      case (IDOMAIN_ELASTIC)
        undef_matpropl(is:ie) = 'elastic'
      end select
      ! default name
      is = (icount-1)*6*MAX_STRING_LEN + (3*MAX_STRING_LEN + 1)
      ie = is + MAX_STRING_LEN
      undef_matpropl(is:ie) = 'tomography_model.xyz'
      ! default tomo-id (unused)
      is = (icount-1)*6*MAX_STRING_LEN + (4*MAX_STRING_LEN + 1)
      ie = is + MAX_STRING_LEN
      write(undef_matpropl(is:ie),*) 0
      ! domain-id
      is = (icount-1)*6*MAX_STRING_LEN + (5*MAX_STRING_LEN + 1)
      ie = is + MAX_STRING_LEN
      write(undef_matpropl(is:ie),*) domain_id
    endif
  enddo
  if (icount /= nundef) stop 'Error icount not equal to ndef'

  ! Boundaries
  call safe_alloc(nodes_ibelm_xmin, ngnod2d, nspec2d_xmin, "nodes_ibelm_xmin")
  call safe_alloc(nodes_ibelm_xmax, ngnod2d, nspec2d_xmax, "nodes_ibelm_xmax")
  call safe_alloc(nodes_ibelm_ymin, ngnod2d, nspec2d_ymin, "nodes_ibelm_ymin")
  call safe_alloc(nodes_ibelm_ymax, ngnod2d, nspec2d_ymax, "nodes_ibelm_ymax")
  call safe_alloc(nodes_ibelm_bottom, ngnod2d, nspec2d_bottom, "nodes_ibelm_bottom")
  call safe_alloc(nodes_ibelm_top, ngnod2d, nspec2d_top, "nodes_ibelm_top")
  call safe_alloc(elmnts_mesh, NGNOD, nspec, "elmnts_mesh")

  do ispec = 1, nspec
    elmnts_mesh(1,ispec) = ibool(1,1,1,ispec)
    elmnts_mesh(2,ispec) = ibool(2,1,1,ispec)
    elmnts_mesh(3,ispec) = ibool(2,2,1,ispec)
    elmnts_mesh(4,ispec) = ibool(1,2,1,ispec)
    elmnts_mesh(5,ispec) = ibool(1,1,2,ispec)
    elmnts_mesh(6,ispec) = ibool(2,1,2,ispec)
    elmnts_mesh(7,ispec) = ibool(2,2,2,ispec)
    elmnts_mesh(8,ispec) = ibool(1,2,2,ispec)
  enddo

  do i=1,nspec2d_xmin
      nodes_ibelm_xmin(1, i) = ibool(1,1,1,ibelm_xmin(i))
      nodes_ibelm_xmin(2, i) = ibool(1,NGLLY_M,1,ibelm_xmin(i))
      nodes_ibelm_xmin(3, i) = ibool(1,1,NGLLZ_M,ibelm_xmin(i))
      nodes_ibelm_xmin(4, i) = ibool(1,NGLLY_M,NGLLZ_M,ibelm_xmin(i))
  enddo
  do i=1,nspec2D_xmax
    nodes_ibelm_xmax(1, i) = ibool(NGLLX_M,1,1,ibelm_xmax(i))
    nodes_ibelm_xmax(2, i) = ibool(NGLLX_M,NGLLY_M,1,ibelm_xmax(i))
    nodes_ibelm_xmax(3, i) = ibool(NGLLX_M,1,NGLLZ_M,ibelm_xmax(i))
    nodes_ibelm_xmax(4, i) = ibool(NGLLX_M,NGLLY_M,NGLLZ_M, ibelm_xmax(i))
  enddo
  do i=1,nspec2D_ymin
    nodes_ibelm_ymin(1, i) = ibool(1,1,1,ibelm_ymin(i))
    nodes_ibelm_ymin(2, i) = ibool(NGLLX_M,1,1,ibelm_ymin(i))
    nodes_ibelm_ymin(3, i) = ibool(1,1,NGLLZ_M,ibelm_ymin(i))
    nodes_ibelm_ymin(4, i) = ibool(NGLLX_M,1,NGLLZ_M,ibelm_ymin(i))
  enddo
  do i=1,nspec2D_ymax
    nodes_ibelm_ymax(1, i) = ibool(NGLLX_M,NGLLY_M,1,ibelm_ymax(i))
    nodes_ibelm_ymax(2, i) = ibool(1,NGLLY_M,1,ibelm_ymax(i))
    nodes_ibelm_ymax(3, i) = ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_ymax(i))
    nodes_ibelm_ymax(4, i) = ibool(1,NGLLY_M,NGLLZ_M,ibelm_ymax(i))
  enddo
  do i=1,NSPEC2D_BOTTOM
    nodes_ibelm_bottom(1, i) = ibool(1,1,1,ibelm_bottom(i))
    nodes_ibelm_bottom(2, i) = ibool(NGLLX_M,1,1,ibelm_bottom(i))
    nodes_ibelm_bottom(3, i) = ibool(NGLLX_M,NGLLY_M,1,ibelm_bottom(i))
    nodes_ibelm_bottom(4, i) = ibool(1,NGLLY_M,1,ibelm_bottom(i))
  enddo
  do i=1,NSPEC2D_TOP
    nodes_ibelm_top(1, i) = ibool(1,1,NGLLZ_M,ibelm_top(i))
    nodes_ibelm_top(2, i) = ibool(NGLLX_M,1,NGLLZ_M,ibelm_top(i))
    nodes_ibelm_top(3, i) = ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_top(i))
    nodes_ibelm_top(4, i) = ibool(1,NGLLY_M,NGLLZ_M,ibelm_top(i))
  enddo

  ! CPML
  call sum_all_i(nspec_CPML,nspec_CPML_total)
  call synchronize_all()
  call bcast_all_singlei(nspec_CPML_total)
  call synchronize_all()

  ! MPI interfaces
  nb_interfaces = 0
  nspec_interfaces_max = 0
  if (NPROC_XI >= 2 .or. NPROC_ETA >= 2) then
    nb_interfaces = 4
    interfaces(W:N) = .true.
    interfaces(NW:SW) = .false.
    if (iproc_xi == 0) then
       nb_interfaces =  nb_interfaces -1
       interfaces(W) = .false.
    endif
    if (iproc_xi == NPROC_XI-1) then
       nb_interfaces =  nb_interfaces -1
       interfaces(E) = .false.
    endif
    if (iproc_eta == 0) then
       nb_interfaces =  nb_interfaces -1
       interfaces(S) = .false.
    endif
    if (iproc_eta == NPROC_ETA-1) then
       nb_interfaces =  nb_interfaces -1
       interfaces(N) = .false.
    endif

    if ((interfaces(W) .eqv. .true.) .and. (interfaces(N) .eqv. .true.)) then
         interfaces(NW) = .true.
         nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(N) .eqv. .true.) .and. (interfaces(E) .eqv. .true.)) then
         interfaces(NE) = .true.
         nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(E) .eqv. .true.) .and. (interfaces(S) .eqv. .true.)) then
         interfaces(SE) = .true.
         nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(W) .eqv. .true.) .and. (interfaces(S) .eqv. .true.)) then
         interfaces(SW) = .true.
         nb_interfaces =  nb_interfaces +1
    endif

    nspec_interface(:) = 0
    if (interfaces(W)) &
        nspec_interface(W) = count(iMPIcut_xi(1,:) .eqv. .true.)
    if (interfaces(E)) &
        nspec_interface(E) = count(iMPIcut_xi(2,:) .eqv. .true.)
    if (interfaces(S)) &
        nspec_interface(S) = count(iMPIcut_eta(1,:) .eqv. .true.)
    if (interfaces(N)) &
        nspec_interface(N) = count(iMPIcut_eta(2,:) .eqv. .true.)
    if (interfaces(NW)) &
        nspec_interface(NW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
    if (interfaces(NE)) &
        nspec_interface(NE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
    if (interfaces(SE)) &
        nspec_interface(SE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))
    if (interfaces(SW)) &
        nspec_interface(SW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))

    nspec_interfaces_max = maxval(nspec_interface)

    call safe_alloc(neighbors_mesh, nb_interfaces, "neighbors_mesh")
    call safe_alloc(num_elmnts_mesh, nb_interfaces, "num_elmnts_mesh")
    call safe_alloc(interfaces_mesh, 6, nspec_interfaces_max, nb_interfaces, "interfaces_mesh")

    interface_num = 1
    neighbors_mesh(:) = 0
    num_elmnts_mesh(:) = 0
    interfaces_mesh(:,:,:) = 0

    if (interfaces(W)) then
      neighbors_mesh(interface_num) = addressing(iproc_xi-1,iproc_eta)
      num_elmnts_mesh(interface_num) = nspec_interface(W)
      ispec_interface = 1
      do ispec = 1,nspec
        if (iMPIcut_xi(1,ispec)) then
          interfaces_mesh(1, ispec_interface, interface_num) = ispec
          interfaces_mesh(2, ispec_interface, interface_num) = 4
          interfaces_mesh(3, ispec_interface, interface_num) &
              = ibool(1,1,1,ispec)
          interfaces_mesh(4, ispec_interface, interface_num) &
              = ibool(1,2,1,ispec)
          interfaces_mesh(5, ispec_interface, interface_num) &
              = ibool(1,1,2,ispec)
          interfaces_mesh(6, ispec_interface, interface_num) &
              = ibool(1,2,2,ispec)
          ispec_interface = ispec_interface + 1
         endif
      enddo
      interface_num = interface_num +1
    endif

    if (interfaces(E)) then
      neighbors_mesh(interface_num) = addressing(iproc_xi+1,iproc_eta)
      num_elmnts_mesh(interface_num) = nspec_interface(E)
      ispec_interface = 1
      do ispec = 1,nspec
        if (iMPIcut_xi(2,ispec)) then
          interfaces_mesh(1, ispec_interface, interface_num) = ispec
          interfaces_mesh(2, ispec_interface, interface_num) = 4
          interfaces_mesh(3, ispec_interface, interface_num) &
              = ibool(2,1,1,ispec)
          interfaces_mesh(4, ispec_interface, interface_num) &
              = ibool(2,2,1,ispec)
          interfaces_mesh(5, ispec_interface, interface_num) &
              = ibool(2,1,2,ispec)
          interfaces_mesh(6, ispec_interface, interface_num) &
              = ibool(2,2,2,ispec)
          ispec_interface = ispec_interface + 1
         endif
      enddo
      interface_num = interface_num +1
    endif

    if (interfaces(S)) then
      neighbors_mesh(interface_num) = addressing(iproc_xi,iproc_eta-1)
      num_elmnts_mesh(interface_num) = nspec_interface(S)
      ispec_interface = 1
      do ispec = 1,nspec
        if (iMPIcut_eta(1,ispec)) then
          interfaces_mesh(1, ispec_interface, interface_num) = ispec
          interfaces_mesh(2, ispec_interface, interface_num) = 4
          interfaces_mesh(3, ispec_interface, interface_num) &
              = ibool(1,1,1,ispec)
          interfaces_mesh(4, ispec_interface, interface_num) &
              = ibool(2,1,1,ispec)
          interfaces_mesh(5, ispec_interface, interface_num) &
              = ibool(1,1,2,ispec)
          interfaces_mesh(6, ispec_interface, interface_num) &
              = ibool(2,1,2,ispec)
          ispec_interface = ispec_interface + 1
        endif
      enddo
      interface_num = interface_num +1
    endif

    if (interfaces(N)) then
      neighbors_mesh(interface_num) = addressing(iproc_xi,iproc_eta+1)
      num_elmnts_mesh(interface_num) = nspec_interface(N)
      ispec_interface = 1
      do ispec = 1,nspec
        if (iMPIcut_eta(2,ispec)) then
          interfaces_mesh(1, ispec_interface, interface_num) = ispec
          interfaces_mesh(2, ispec_interface, interface_num) = 4
          interfaces_mesh(3, ispec_interface, interface_num) &
              = ibool(2,2,1,ispec)
          interfaces_mesh(4, ispec_interface, interface_num) &
              = ibool(1,2,1,ispec)
          interfaces_mesh(5, ispec_interface, interface_num) &
              = ibool(2,2,2,ispec)
          interfaces_mesh(6, ispec_interface, interface_num) &
              = ibool(1,2,2,ispec)
          ispec_interface = ispec_interface + 1
        endif
      enddo
      interface_num = interface_num +1
    endif

    if (interfaces(NW)) then
      neighbors_mesh(interface_num) = addressing(iproc_xi-1,iproc_eta+1)
      num_elmnts_mesh(interface_num) = nspec_interface(NW)
      ispec_interface = 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          interfaces_mesh(1, ispec_interface, interface_num) = ispec
          interfaces_mesh(2, ispec_interface, interface_num) = 2
          interfaces_mesh(3, ispec_interface, interface_num) &
              = ibool(1,2,1,ispec)
          interfaces_mesh(4, ispec_interface, interface_num) &
              = ibool(1,2,2,ispec)
          interfaces_mesh(5, ispec_interface, interface_num) = -1
          interfaces_mesh(6, ispec_interface, interface_num) = -1
          ispec_interface = ispec_interface + 1
        endif
      enddo
      interface_num = interface_num +1
    endif

    if (interfaces(NE)) then
      neighbors_mesh(interface_num) = addressing(iproc_xi+1,iproc_eta+1)
      num_elmnts_mesh(interface_num) = nspec_interface(NE)
      ispec_interface = 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          interfaces_mesh(1, ispec_interface, interface_num) = ispec
          interfaces_mesh(2, ispec_interface, interface_num) = 2
          interfaces_mesh(3, ispec_interface, interface_num) &
              = ibool(2,2,1,ispec)
          interfaces_mesh(4, ispec_interface, interface_num) &
              = ibool(2,2,2,ispec)
          interfaces_mesh(5, ispec_interface, interface_num) = -1
          interfaces_mesh(6, ispec_interface, interface_num) = -1
          ispec_interface = ispec_interface + 1
        endif
      enddo
      interface_num = interface_num +1
    endif

    if (interfaces(SE)) then
      neighbors_mesh(interface_num) = addressing(iproc_xi+1,iproc_eta-1)
      num_elmnts_mesh(interface_num) = nspec_interface(SE)
      ispec_interface = 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          interfaces_mesh(1, ispec_interface, interface_num) = ispec
          interfaces_mesh(2, ispec_interface, interface_num) = 2
          interfaces_mesh(3, ispec_interface, interface_num) &
              = ibool(2,1,1,ispec)
          interfaces_mesh(4, ispec_interface, interface_num) &
              = ibool(2,1,2,ispec)
          interfaces_mesh(5, ispec_interface, interface_num) = -1
          interfaces_mesh(6, ispec_interface, interface_num) = -1
          ispec_interface = ispec_interface + 1
        endif
      enddo
      interface_num = interface_num +1
    endif

    if (interfaces(SW)) then
      neighbors_mesh(interface_num) = addressing(iproc_xi-1,iproc_eta-1)
      num_elmnts_mesh(interface_num) = nspec_interface(SW)
      ispec_interface = 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          interfaces_mesh(1, ispec_interface, interface_num) = ispec
          interfaces_mesh(2, ispec_interface, interface_num) = 2
          interfaces_mesh(3, ispec_interface, interface_num) &
              = ibool(1,1,1,ispec)
          interfaces_mesh(4, ispec_interface, interface_num) &
              = ibool(1,1,2,ispec)
          interfaces_mesh(5, ispec_interface, interface_num) = -1
          interfaces_mesh(6, ispec_interface, interface_num) = -1
          ispec_interface = ispec_interface + 1
        endif
      enddo
      interface_num = interface_num +1
    endif
  else
    ! only one MPI process
    call safe_alloc(neighbors_mesh, nb_interfaces, "neighbors_mesh")
    call safe_alloc(num_elmnts_mesh, nb_interfaces, "num_elmnts_mesh")
    call safe_alloc(interfaces_mesh, 6, nspec_interfaces_max, nb_interfaces, "interfaces_mesh")
  endif

  !-----------------------------------------------------------------.
  ! Get maximum value for each variable used to define a local_dim. |
  ! ADIOS write equally sized chunks for each processor.            |
  !-----------------------------------------------------------------'

  ! Filling a temporary array to avoid doing allreduces for each var.
  max_global_values(1) = nglob
  max_global_values(2) = nspec
  max_global_values(3) = NMATERIALS
  max_global_values(4) = nspec2d_xmin
  max_global_values(5) = nspec2d_xmax
  max_global_values(6) = nspec2d_ymin
  max_global_values(7) = nspec2d_ymax
  max_global_values(8) = nspec2d_bottom
  max_global_values(9) = nspec2d_top
  max_global_values(10) = nb_interfaces
  max_global_values(11) = nspec_interfaces_max

  call max_allreduce_i(max_global_values,num_vars)

  nglob_wmax          = max_global_values(1)
  nspec_wmax          = max_global_values(2)
  nmaterials_wmax     = max_global_values(3)
  nspec2d_xmin_wmax   = max_global_values(4)
  nspec2d_xmax_wmax   = max_global_values(5)
  nspec2d_ymin_wmax   = max_global_values(6)
  nspec2d_ymax_wmax   = max_global_values(7)
  nspec2d_bottom_wmax = max_global_values(8)
  nspec2d_top_wmax    = max_global_values(9)
  nb_interfaces_wmax  = max_global_values(10)
  nspec_interfaces_max_wmax = max_global_values(11)

  !-----------------------------------.
  ! Setup ADIOS for the current group |
  !-----------------------------------'
  groupsize = 0
  output_name = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // "/Database.bp"
  call adios_declare_group(group, group_name, '', 1, ier)
  call adios_select_method(group, ADIOS_TRANSPORT_METHOD, '', '', ier)

  !------------------------.
  ! Define ADIOS Variables |
  !------------------------'
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(ngllx))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nglly))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(ngllz))

  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(ngllx_m))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nglly_m))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(ngllz_m))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(ngnod))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(ngnod2d))

  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nglob))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nspec))

  call define_adios_scalar(group, groupsize, '', "nmaterials",ndef)
  call define_adios_scalar(group, groupsize, '', "nundef_materials",nundef)

  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nspec2d_xmin))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nspec2d_xmax))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nspec2d_ymin))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nspec2d_ymax))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nspec2d_bottom))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nspec2d_top))

  call define_adios_scalar(group, groupsize, '', "nspec_cpml_total",nspec_cpml_total)
  if (nspec_cpml_total > 0) call define_adios_scalar(group, groupsize, '', "nspec_cpml",nspec_cpml)

  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nb_interfaces))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nspec_interfaces_max))

  local_dim = 3 * nglob_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(nodes_coords))

  if (ndef /= 0) then
    local_dim = 16 * ndef
    call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(matpropl))
  endif

  if (nundef /= 0) then
    local_dim = 6 * nundef * MAX_STRING_LEN
    ! all processes will have an identical entry, uses "local" variable to store string
    call define_adios_local_string_1D_array(group, groupsize, local_dim, '', "undef_mat_prop", undef_matpropl(1:local_dim))
    ! note: global arrays with strings doesn't work yet... (segmentation fault)
    !call define_adios_global_array1D(group, groupsize, local_dim, '', "undef_mat_prop",undef_matpropl)
  endif

  local_dim = 2 * nspec_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(material_index))

  local_dim = NGLLX_M * NGLLY_M * NGLLZ_M * nspec_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(elmnts_mesh))

  local_dim = nspec2d_xmin_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(ibelm_xmin))
  local_dim = nspec2d_xmax_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(ibelm_xmax))
  local_dim = nspec2d_ymin_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(ibelm_ymin))
  local_dim = nspec2d_ymax_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(ibelm_ymax))
  local_dim = nspec2d_bottom_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(ibelm_bottom))
  local_dim = nspec2d_top_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(ibelm_top))

  local_dim = NGNOD2D * nspec2d_xmin_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(nodes_ibelm_xmin))
  local_dim = NGNOD2D * nspec2d_xmax_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(nodes_ibelm_xmax))
  local_dim = NGNOD2D * nspec2d_ymin_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(nodes_ibelm_ymin))
  local_dim = NGNOD2D * nspec2d_ymax_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(nodes_ibelm_ymax))
  local_dim = NGNOD2D * nspec2d_bottom_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(nodes_ibelm_bottom))
  local_dim = NGNOD2D * nspec2d_top_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(nodes_ibelm_top))

  if (nspec_CPML_total > 0) then
     local_dim=nspec_CPML
     call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(CPML_to_spec))
     call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(CPML_regions))
     local_dim=nspec
     call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(is_CPML))
  endif

  if (nb_interfaces_wmax > 0) then
    local_dim = nb_interfaces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(neighbors_mesh))
    call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(num_elmnts_mesh))
    local_dim = 6 * nb_interfaces_wmax * nspec_interfaces_max_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, '', STRINGIFY_VAR(interfaces_mesh))
  endif

  !------------------------------------------------------------.
  ! Open an handler to the ADIOS file and setup the group size |
  !------------------------------------------------------------'
  call world_get_comm(comm)

  call adios_open(handle, group_name, output_name, "w", comm, ier);
  call adios_group_size (handle, groupsize, totalsize, ier)

  !------------------------------------------.
  ! Write previously defined ADIOS variables |
  !------------------------------------------'
  call adios_write(handle, STRINGIFY_VAR(ngllx), ier)
  call adios_write(handle, STRINGIFY_VAR(nglly), ier)
  call adios_write(handle, STRINGIFY_VAR(ngllz), ier)

  call adios_write(handle, STRINGIFY_VAR(ngllx_m), ier)
  call adios_write(handle, STRINGIFY_VAR(nglly_m), ier)
  call adios_write(handle, STRINGIFY_VAR(ngllz_m), ier)
  call adios_write(handle, STRINGIFY_VAR(ngnod), ier)
  call adios_write(handle, STRINGIFY_VAR(ngnod2d), ier)

  call adios_write(handle, STRINGIFY_VAR(nglob), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec), ier)

  call adios_write(handle, "nmaterials",ndef, ier)
  call adios_write(handle, "nundef_materials",nundef, ier)

  call adios_write(handle, STRINGIFY_VAR(nspec2d_xmin), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec2d_xmax), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec2d_ymin), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec2d_ymax), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec2d_bottom), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec2d_top), ier)

  call adios_write(handle, "nspec_cpml_total", nspec_cpml_total, ier)
  if (nspec_CPML_total > 0) call adios_write(handle, "nspec_cpml", nspec_cpml, ier)

  call adios_write(handle, STRINGIFY_VAR(nb_interfaces), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec_interfaces_max), ier)

  ! NOTE: Do not put any wmax variables, it will try to access
  !       too many values in the arrays.
  local_dim = 3 * nglob
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   "nodes_coords", transpose(nodes_coords))

  ! materials
  if (ndef /= 0) then
    local_dim = 16 * ndef
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(matpropl))
  endif
  if (nundef /= 0) then
    local_dim = 6 * nundef * MAX_STRING_LEN
    ! "local" variable, all process will write the same entry
    call adios_write(handle, "undef_mat_prop", undef_matpropl(1:local_dim),ier)

    ! note: global arrays don't work yet with strings, produces segmentation faults (ADIOS error?)...
    !call write_adios_global_string_1d_array(handle, myrank, sizeprocs, local_dim, local_dim*sizeprocs, local_dim*myrank, &
    !                                        "undef_mat_prop",undef_matpropl(1:local_dim))
  endif

  local_dim = 2 * nspec
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(material_index))
  ! WARNING: the order is a little bit different than for Fortran output
  !          It should not matter, but it may.
  local_dim = NGLLX_M * NGLLY_M * NGLLZ_M * nspec_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(elmnts_mesh))

  ! model boundary surfaces
  local_dim = nspec2d_xmin_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ibelm_xmin))
  local_dim = nspec2d_xmax_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ibelm_xmax))
  local_dim = nspec2d_ymin_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ibelm_ymin))
  local_dim = nspec2d_ymax_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ibelm_ymax))
  local_dim = nspec2d_bottom_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ibelm_bottom))
  local_dim = nspec2d_top_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ibelm_top))

  ! needs surface boundary points
  if (nspec2d_xmin > 0) then
    local_dim = NGNOD2D * nspec2d_xmin_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(nodes_ibelm_xmin))
  endif

  if (nspec2d_xmax > 0) then
    local_dim = NGNOD2D * nspec2d_xmax_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(nodes_ibelm_xmax))
  endif

  if (nspec2d_ymin > 0) then
    local_dim = NGNOD2D * nspec2d_ymin_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(nodes_ibelm_ymin))
  endif

  if (nspec2d_ymax > 0) then
    local_dim = NGNOD2D * nspec2d_ymax_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(nodes_ibelm_ymax))
  endif

  if (nspec2d_bottom > 0) then
    local_dim = NGNOD2D * nspec2d_bottom_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(nodes_ibelm_bottom))
  endif

  if (nspec2d_top > 0) then
    local_dim = NGNOD2D * nspec2d_top_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(nodes_ibelm_top))
  endif

  ! CPML
  if (nspec_CPML_total > 0) then
     local_dim=nspec_CPML
     call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(CPML_to_spec))
     call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(CPML_regions))

     local_dim=nspec
     call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(is_CPML))
  endif

  ! MPI interfaces
  if (nb_interfaces > 0) then
    local_dim = nb_interfaces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(neighbors_mesh))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(num_elmnts_mesh))
    local_dim = 6 * nb_interfaces_wmax * nspec_interfaces_max_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(interfaces_mesh))
  endif

  !----------------------------------.
  ! Perform the actual write to disk |
  !----------------------------------'
  call adios_set_path(handle, '', ier)
  call adios_close(handle, ier)

  !---------------------------.
  ! Clean up temporary arrays |
  !---------------------------'
  call safe_dealloc(nodes_ibelm_xmin, "nodes_ibelm_xmin")
  call safe_dealloc(nodes_ibelm_xmax, "nodes_ibelm_xmax")
  call safe_dealloc(nodes_ibelm_ymin, "nodes_ibelm_ymin")
  call safe_dealloc(nodes_ibelm_ymax, "nodes_ibelm_ymax")
  call safe_dealloc(nodes_ibelm_bottom, "nodes_ibelm_bottom")
  call safe_dealloc(nodes_ibelm_top, "nodes_ibelm_top")

  call safe_dealloc(neighbors_mesh, "neighbors_mesh")
  call safe_dealloc(num_elmnts_mesh, "num_elmnts_mesh")
  call safe_dealloc(interfaces_mesh, "interfaces_mesh")
  call safe_dealloc(elmnts_mesh, "elmnts_mesh")

  call safe_dealloc(matpropl,"matpropl")

end subroutine save_databases_adios
