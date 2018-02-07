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

  subroutine save_databases(prname,nspec,nglob,iproc_xi,iproc_eta, &
                            NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta, &
                            ibool,nodes_coords,ispec_material_id, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                            NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            NMATERIALS,material_properties, &
                            nspec_CPML,CPML_to_spec,CPML_regions,is_CPML, &
                            xstore, ystore, zstore)

  use constants, only: MAX_STRING_LEN,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,SAVE_MESH_AS_CUBIT
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE

  implicit none

  include "constants_meshfem3D.h"

  ! number of spectral elements in each block
  integer nspec

  ! number of vertices in each block
  integer nglob

  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8
  integer iproc_xi,iproc_eta
  integer NPROC_XI,NPROC_ETA
  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)
  integer addressing(0:NPROC_XI-1,0:NPROC_ETA-1)

  ! arrays with the mesh
  integer ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)
  double precision :: nodes_coords(nglob,3)

  !! VM VM add all GLL points for Axisem coupling
  double precision, dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec) :: xstore, ystore, zstore

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
  ! second dimension : #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id #material_id
  double precision , dimension(NMATERIALS,7) ::  material_properties

  ! CPML
  integer, intent(in) :: nspec_CPML
  integer, dimension(nspec_CPML), intent(in) :: CPML_to_spec,CPML_regions
  logical, dimension(nspec), intent(in) :: is_CPML
  integer :: nspec_CPML_total,ispec_CPML

  double precision , dimension(16) :: matpropl
  integer :: i,ispec,iglob,ier

  ! name of the database files
  character(len=MAX_STRING_LEN) :: prname

  ! for MPI interfaces
  integer ::  nb_interfaces,nspec_interfaces_max
  logical, dimension(8) ::  interfaces
  integer, dimension(8) ::  nspec_interface

  integer, parameter :: IIN_database = 15

  integer :: ndef,nundef
  integer :: mat_id,domain_id
  ! there was a potential bug here if nspec is big
  !integer,dimension(2,nspec) :: material_index
  integer,dimension(:,:),allocatable :: material_index
  character(len=MAX_STRING_LEN), dimension(6,1) :: undef_mat_prop

  ! assignes material index
  ! format: (1,ispec) = #material_id , (2,ispec) = #material_definition
  allocate(material_index(2,nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating array material_index'
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
    ! material properties format: #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id #material_id
    mat_id = material_properties(i,7)
    if (mat_id > 0) ndef = ndef + 1
    if (mat_id < 0) nundef = nundef + 1
  enddo
  !debug
  !print *,'materials def/undef: ',ndef,nundef

  ! opens database file
  open(unit=IIN_database,file=prname(1:len_trim(prname))//'Database', &
        status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening Database file: ',prname(1:len_trim(prname))//'Database'
    stop 'error opening Database file'
  endif

  ! global nodes
  write(IIN_database) nglob
  do iglob=1,nglob
     write(IIN_database) iglob,nodes_coords(iglob,1),nodes_coords(iglob,2),nodes_coords(iglob,3)
  enddo

  ! materials
  ! format: #number of defined materials #number of undefined materials
  write(IIN_database) ndef, nundef

  ! writes out defined materials
  do i = 1,NMATERIALS
    ! material properties format: #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id #material_id
    mat_id = material_properties(i,7)
    if (mat_id > 0) then
      ! pad dummy zeros to fill up 16 entries (poroelastic medium not allowed)
      matpropl(:) = 0.d0
      ! material properties format: #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id
      matpropl(1:6) = material_properties(i,1:6)
      ! fills adds arbitrary value for Q_kappa
      matpropl(7) = 9999.0
      write(IIN_database) matpropl(:)
    endif
  enddo

  ! writes out undefined materials
  do i = 1,NMATERIALS
    domain_id = material_properties(i,6)
    mat_id = material_properties(i,7)
    if (mat_id < 0) then
      ! format:
      ! #material_id #type-keyword #domain-name #tomo-filename #tomo_id #domain_id
      ! format example tomography: -1 tomography elastic tomography_model.xyz 0 2
      undef_mat_prop(:,:) = ''
      ! material id
      write(undef_mat_prop(1,1),*) mat_id
      ! name
      undef_mat_prop(2,1) = 'tomography'
      select case (domain_id)
      case (IDOMAIN_ACOUSTIC)
        undef_mat_prop(3,1) = 'acoustic'
      case (IDOMAIN_ELASTIC)
        undef_mat_prop(3,1) = 'elastic'
      end select
      ! default name
      undef_mat_prop(4,1) = 'tomography_model.xyz'
      ! default tomo-id (unused)
      write(undef_mat_prop(5,1),*) 0
      ! domain-id
      write(undef_mat_prop(6,1),*) domain_id
      ! debug
      !print *,'undef mat: ',undef_mat_prop
      ! writes out properties
      write(IIN_database) undef_mat_prop
    endif
  enddo

  ! spectral elements
  write(IIN_database) nspec
  do ispec=1,nspec
      write(IIN_database) ispec,material_index(1,ispec),material_index(2,ispec), &
           ibool(1,1,1,ispec),ibool(2,1,1,ispec),ibool(2,2,1,ispec),ibool(1,2,1,ispec),ibool(1,1,2,ispec), &
           ibool(2,1,2,ispec),ibool(2,2,2,ispec),ibool(1,2,2,ispec)
  enddo

  ! Boundaries
  write(IIN_database) 1,nspec2D_xmin
  write(IIN_database) 2,nspec2D_xmax
  write(IIN_database) 3,nspec2D_ymin
  write(IIN_database) 4,nspec2D_ymax
  write(IIN_database) 5,NSPEC2D_BOTTOM
  write(IIN_database) 6,NSPEC2D_TOP

  do i=1,nspec2D_xmin
     write(IIN_database) ibelm_xmin(i),ibool(1,1,1,ibelm_xmin(i)),ibool(1,NGLLY_M,1,ibelm_xmin(i)), &
          ibool(1,1,NGLLZ_M,ibelm_xmin(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_xmin(i))
  enddo
  do i=1,nspec2D_xmax
     write(IIN_database) ibelm_xmax(i),ibool(NGLLX_M,1,1,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,1,ibelm_xmax(i)), &
          ibool(NGLLX_M,1,NGLLZ_M,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_xmax(i))
  enddo
  do i=1,nspec2D_ymin
     write(IIN_database) ibelm_ymin(i),ibool(1,1,1,ibelm_ymin(i)),ibool(NGLLX_M,1,1,ibelm_ymin(i)), &
          ibool(1,1,NGLLZ_M,ibelm_ymin(i)),ibool(NGLLX_M,1,NGLLZ_M,ibelm_ymin(i))
  enddo
  do i=1,nspec2D_ymax
     write(IIN_database) ibelm_ymax(i),ibool(NGLLX_M,NGLLY_M,1,ibelm_ymax(i)),ibool(1,NGLLY_M,1,ibelm_ymax(i)), &
          ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_ymax(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_ymax(i))
  enddo
  do i=1,NSPEC2D_BOTTOM
     write(IIN_database) ibelm_bottom(i),ibool(1,1,1,ibelm_bottom(i)),ibool(NGLLX_M,1,1,ibelm_bottom(i)), &
          ibool(NGLLX_M,NGLLY_M,1,ibelm_bottom(i)),ibool(1,NGLLY_M,1,ibelm_bottom(i))
  enddo
  do i=1,NSPEC2D_TOP
     write(IIN_database) ibelm_top(i),ibool(1,1,NGLLZ_M,ibelm_top(i)),ibool(NGLLX_M,1,NGLLZ_M,ibelm_top(i)), &
          ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_top(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_top(i))
  enddo

  ! CPML
  call sum_all_i(nspec_CPML,nspec_CPML_total)
  call synchronize_all()
  call bcast_all_singlei(nspec_CPML_total)
  call synchronize_all()

  write(IIN_database) nspec_CPML_total
  if (nspec_CPML_total > 0) then
     write(IIN_database) nspec_CPML
     do ispec_CPML=1,nspec_CPML
        write(IIN_database) CPML_to_spec(ispec_CPML), CPML_regions(ispec_CPML)
     enddo
     do ispec=1,nspec
        write(IIN_database) is_CPML(ispec)
     enddo
  endif

  ! MPI Interfaces
  if (NPROC_XI >= 2 .or. NPROC_ETA >= 2) then
    ! determines number of MPI interfaces for each slice
    nb_interfaces = 4
    interfaces(W:N) = .true.
    interfaces(NW:SW) = .false.

    ! slices at model boundaries
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

    ! slices in middle of model
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
    if (interfaces(W))  nspec_interface(W) = count(iMPIcut_xi(1,:) .eqv. .true.)
    if (interfaces(E))  nspec_interface(E) = count(iMPIcut_xi(2,:) .eqv. .true.)
    if (interfaces(S))  nspec_interface(S) = count(iMPIcut_eta(1,:) .eqv. .true.)
    if (interfaces(N))  nspec_interface(N) = count(iMPIcut_eta(2,:) .eqv. .true.)
    if (interfaces(NW))  nspec_interface(NW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
    if (interfaces(NE))  nspec_interface(NE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
    if (interfaces(SE))  nspec_interface(SE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))
    if (interfaces(SW))  nspec_interface(SW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))

    nspec_interfaces_max = maxval(nspec_interface)

    write(IIN_database) nb_interfaces,nspec_interfaces_max

    if (interfaces(W)) then
      write(IIN_database) addressing(iproc_xi-1,iproc_eta),nspec_interface(W)
      do ispec = 1,nspec
        if (iMPIcut_xi(1,ispec)) write(IIN_database) ispec,4,ibool(1,1,1,ispec),ibool(1,2,1,ispec), &
                                                    ibool(1,1,2,ispec),ibool(1,2,2,ispec)
      enddo
    endif

    if (interfaces(E)) then
      write(IIN_database) addressing(iproc_xi+1,iproc_eta),nspec_interface(E)
      do ispec = 1,nspec
        if (iMPIcut_xi(2,ispec)) write(IIN_database) ispec,4,ibool(2,1,1,ispec),ibool(2,2,1,ispec), &
                                                    ibool(2,1,2,ispec),ibool(2,2,2,ispec)
      enddo
    endif

    if (interfaces(S)) then
      write(IIN_database) addressing(iproc_xi,iproc_eta-1),nspec_interface(S)
      do ispec = 1,nspec
        if (iMPIcut_eta(1,ispec)) write(IIN_database) ispec,4,ibool(1,1,1,ispec),ibool(2,1,1,ispec), &
                                                     ibool(1,1,2,ispec),ibool(2,1,2,ispec)
      enddo
    endif

    if (interfaces(N)) then
      write(IIN_database) addressing(iproc_xi,iproc_eta+1),nspec_interface(N)
      do ispec = 1,nspec
        if (iMPIcut_eta(2,ispec)) write(IIN_database) ispec,4,ibool(2,2,1,ispec),ibool(1,2,1,ispec), &
                                                     ibool(2,2,2,ispec),ibool(1,2,2,ispec)
      enddo
    endif

    if (interfaces(NW)) then
      write(IIN_database) addressing(iproc_xi-1,iproc_eta+1),nspec_interface(NW)
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          write(IIN_database) ispec,2,ibool(1,2,1,ispec),ibool(1,2,2,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(NE)) then
      write(IIN_database) addressing(iproc_xi+1,iproc_eta+1),nspec_interface(NE)
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          write(IIN_database) ispec,2,ibool(2,2,1,ispec),ibool(2,2,2,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(SE)) then
      write(IIN_database) addressing(iproc_xi+1,iproc_eta-1),nspec_interface(SE)
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          write(IIN_database) ispec,2,ibool(2,1,1,ispec),ibool(2,1,2,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(SW)) then
      write(IIN_database) addressing(iproc_xi-1,iproc_eta-1),nspec_interface(SW)
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          write(IIN_database) ispec,2,ibool(1,1,1,ispec),ibool(1,1,2,ispec),-1,-1
        endif
      enddo
    endif

  else

    ! only one slice, no MPI interfaces
    write(IIN_database) 0,0

  endif

  close(IIN_database)

  ! CUBIT output
  if (SAVE_MESH_AS_CUBIT) then
    ! only for single process at the moment
    if (NPROC_XI == 1 .and. NPROC_ETA == 1) then
      !! VM VM add outputs as CUBIT
      call save_output_mesh_files_as_cubit(nspec,nglob, ibool,nodes_coords, ispec_material_id, &
                                           nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NMATERIALS,material_properties, &
                                           ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top)
      ! output for AxiSEM coupling
      if (COUPLE_WITH_INJECTION_TECHNIQUE) then
        call save_output_mesh_files_for_coupled_model(nspec, &
                                           nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                                           ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                           xstore,ystore,zstore)
      endif
    endif
  endif

  deallocate(material_index)

  end subroutine save_databases

!---------------------------------------------------------------------------------------------------------------

  !! VM VM add subroutine to save meshes in case of a single MPI process
  subroutine save_output_mesh_files_as_cubit(nspec,nglob, ibool,nodes_coords, ispec_material_id, &
                                     nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                     NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NMATERIALS,material_properties, &
                                     ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top)

  use constants, only: MAX_STRING_LEN,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC, NGLLX, NGLLY, NGLLZ, NDIM, ZERO

  implicit none

  include "constants_meshfem3D.h"

  integer, parameter :: IIN_database = 15
  ! number of spectral elements in each block
  integer :: nspec

  ! number of vertices in each block
  integer :: nglob

  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8

  ! arrays with the mesh
  integer :: ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)
  double precision :: nodes_coords(nglob,3)

  integer :: ispec_material_id(nspec)

  ! boundary parameters locator
  integer :: NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer :: ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer :: ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer :: ibelm_bottom(NSPEC2D_BOTTOM)
  integer :: ibelm_top(NSPEC2D_TOP)

  ! material properties
  integer :: NMATERIALS
  ! first dimension  : material_id
  ! second dimension : #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id #material_id
  double precision , dimension(NMATERIALS,7) ::  material_properties

  ! local parameters
  integer :: i,ispec,iglob,ier
  integer :: mat_id,domain_id
  double precision  :: z_bottom

  z_bottom = 0. ! will shift coordinates in z-direction

  ! outputs mesh as files in MESH/ directory
  ! (used for CUBIT models stored in a specfem-readable way; will need to run xdecompose_mesh with these files to continue)

  open(IIN_database, file = 'MESH/nummaterial_velocity_file',status='unknown',action='write',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ','MESH/nummaterial_velocity_file'
    print *,'Please check if directory MESH/ exists for saving mesh files as cubit...'
    stop 'Error opening mesh file'
  endif

  do i = 1,NMATERIALS
     domain_id = material_properties(i,6)
     mat_id =  material_properties(i,7)
     if ( domain_id > 0) then
        write(IIN_database,'(2i6,5f15.5,i6)') domain_id,mat_id,material_properties(i,1:3),9999.,9999.,0
     else
       write(*,*) 'STOP: undefined mat not yet implemented'
       stop
     endif
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/materials_file')
  do ispec = 1, nspec
     write(IIN_database,*) ispec,ispec_material_id(ispec)
  enddo

  open(IIN_database,file='MESH/nodes_coords_file')
  write(IIN_database,*) nglob
  do iglob=1,nglob
     write(IIN_database,'(i14,3x,3(f20.5,1x))') iglob,nodes_coords(iglob,1),nodes_coords(iglob,2), &
          nodes_coords(iglob,3)-z_bottom
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/mesh_file')
  write(IIN_database,*) nspec
  do ispec=1,nspec
     write(IIN_database,'(9i15)')  ispec,ibool(1,1,1,ispec),ibool(2,1,1,ispec), &
          ibool(2,2,1,ispec),ibool(1,2,1,ispec), &
          ibool(1,1,2,ispec),ibool(2,1,2,ispec), &
          ibool(2,2,2,ispec),ibool(1,2,2,ispec)
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_xmin')
  write(IIN_database,*) nspec2D_xmin
  do i=1,nspec2D_xmin
     write(IIN_database,'(5(i10,1x))') ibelm_xmin(i),ibool(1,1,1,ibelm_xmin(i)),ibool(1,NGLLY_M,1,ibelm_xmin(i)), &
        ibool(1,1,NGLLZ_M,ibelm_xmin(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_xmin(i))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_xmax')
  write(IIN_database,*) nspec2D_xmax
  do i=1,nspec2D_xmax
     write(IIN_database,'(5(i10,1x))') ibelm_xmax(i),ibool(NGLLX_M,1,1,ibelm_xmax(i)), &
          ibool(NGLLX_M,NGLLY_M,1,ibelm_xmax(i)), ibool(NGLLX_M,1,NGLLZ_M,ibelm_xmax(i)), &
          ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_xmax(i))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_ymin')
  write(IIN_database,*) nspec2D_ymin
  do i=1,nspec2D_ymin
     write(IIN_database,'(5(i10,1x))') ibelm_ymin(i),ibool(1,1,1,ibelm_ymin(i)),ibool(NGLLX_M,1,1,ibelm_ymin(i)), &
          ibool(1,1,NGLLZ_M,ibelm_ymin(i)),ibool(NGLLX_M,1,NGLLZ_M,ibelm_ymin(i))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_ymax')
  write(IIN_database,*) nspec2D_ymax
  do i=1,nspec2D_ymax
     write(IIN_database,'(5(i10,1x))') ibelm_ymax(i),ibool(NGLLX_M,NGLLY_M,1,ibelm_ymax(i)), &
          ibool(1,NGLLY_M,1,ibelm_ymax(i)), ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_ymax(i)), &
          ibool(1,NGLLY_M,NGLLZ_M,ibelm_ymax(i))
  enddo


  open(IIN_database,file='MESH/absorbing_surface_file_bottom')
  write(IIN_database,*) NSPEC2D_BOTTOM
  do i=1,NSPEC2D_BOTTOM
     write(IIN_database,'(5(i10,1x))') ibelm_bottom(i),ibool(1,1,1,ibelm_bottom(i)),ibool(NGLLX_M,1,1,ibelm_bottom(i)), &
          ibool(NGLLX_M,NGLLY_M,1,ibelm_bottom(i)),ibool(1,NGLLY_M,1,ibelm_bottom(i))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/free_or_absorbing_surface_file_zmax')
  write(IIN_database,*) NSPEC2D_TOP
  do i=1,NSPEC2D_TOP
     write(IIN_database,'(5(i10,1x))') ibelm_top(i),ibool(1,1,NGLLZ_M,ibelm_top(i)), &
          ibool(NGLLX_M,1,NGLLZ_M,ibelm_top(i)),ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_top(i)), &
          ibool(1,NGLLY_M,NGLLZ_M,ibelm_top(i))
  enddo
  close(IIN_database)

  end subroutine save_output_mesh_files_as_cubit



!---------------------------------------------------------------------------------------------------------------

  !! VM VM add subroutine to save meshes in case of a single MPI process
  subroutine save_output_mesh_files_for_coupled_model(nspec, &
                                              nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                                              ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                              xgrid,ygrid,zgrid)

  use constants, only: MAX_STRING_LEN,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC, NGLLX, NGLLY, NGLLZ, NDIM, ZERO, &
    INJECTION_TECHNIQUE_IS_AXISEM

  use shared_parameters, only: NGNOD,COUPLE_WITH_INJECTION_TECHNIQUE,INJECTION_TECHNIQUE_TYPE

  implicit none

  include "constants_meshfem3D.h"

  ! number of spectral elements in each block
  integer :: nspec

  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8

  !! VM VM add all GLL points for Axisem coupling
  double precision, dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec) :: xgrid, ygrid, zgrid

  ! boundary parameters locator
  integer :: NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer :: ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer :: ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer :: ibelm_bottom(NSPEC2D_BOTTOM)
  integer :: ibelm_top(NSPEC2D_TOP)

  integer :: i,ispec

  double precision  :: z_bottom

  ! for axisem coupling case  ( only serial case for mesher use scotch after)
  integer, parameter :: myrank = 0
  integer, parameter :: nlayer = 12 !! (number of layer in the model iasp91, or ak135, or prem (one more layer than the model)
  double precision, parameter :: GAUSSALPHA = 0.d0, GAUSSBETA = 0.d0
  double precision   :: rotation_matrix(3,3)
  double precision   :: zlayer(nlayer), vpv(nlayer,4), vsv(nlayer,4), density(nlayer,4)
  integer            :: ilayer, updown(NGLLZ)

  !! GLL points
  double precision, dimension(:,:,:), allocatable  ::  longitud, latitud, radius
  double precision, dimension(:,:,:), allocatable  ::  xstore, ystore, zstore
   !! Element control points
  double precision, dimension(:), allocatable :: xelm, yelm, zelm

  !! 3D shape functions and their derivatives
  double precision, dimension(:,:,:,:), allocatable    :: shape3D
  double precision, dimension(:,:,:,:,:), allocatable  :: dershape3D
  !! GLL points and weights of integration
  double precision, dimension(:), allocatable  :: xigll, yigll, zigll, wxgll, wygll, wzgll

  double precision  :: deg2rad
  double precision  :: ANGULAR_WIDTH_ETA_RAD, ANGULAR_WIDTH_XI_RAD
  double precision  :: lat_center_chunk, lon_center_chunk, chunk_depth, chunk_azi
  double precision  :: radius_of_box_top

  integer :: ielm, j,k, imin,imax,jmin,jmax,kmin,kmax
  integer :: nel_lat, nel_lon, nel_depth
  logical :: buried_box

  character(len=10)  :: line
  character(len=250) :: model1D_file

  ! safety check
  if (.not. COUPLE_WITH_INJECTION_TECHNIQUE) return

!! VM VM add files in case of AxiSEM coupling
  if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then

1000 format(3f30.10)

    z_bottom = 0. ! will shift coordinates in z-direction

    allocate(longitud(NGLLX,NGLLY,NGLLZ), latitud(NGLLX,NGLLY,NGLLZ), radius(NGLLX,NGLLY,NGLLZ))
    allocate(xstore(NGLLX,NGLLY,NGLLZ), ystore(NGLLX,NGLLY,NGLLZ), zstore(NGLLX,NGLLY,NGLLZ))
    allocate(xelm(NGNOD), yelm(NGNOD), zelm(NGNOD))
    allocate(xigll(NGLLX), yigll(NGLLY), zigll(NGLLZ), wxgll(NGLLX),wygll(NGLLY), wzgll(NGLLZ))

    deg2rad = 3.141592653589793d0/180.d0

    !
    !--- set up coordinates of the Gauss-Lobatto-Legendre points
    !

    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
    call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

    !
    !--- if number of points is odd, the middle abscissa is exactly zero
    !
    if (mod(NGLLX,2) /= 0) xigll((NGLLX - 1)/2 + 1) = ZERO
    if (mod(NGLLY,2) /= 0) yigll((NGLLY - 1)/2 + 1) = ZERO
    if (mod(NGLLZ,2) /= 0) zigll((NGLLZ - 1)/2 + 1) = ZERO

    !
    !--- get the 3-D shape functions
    !
    allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ),dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ))
    call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll,NGNOD)
    !

    !! reading parameters for coupling

    open(90, file='MESH/ParFileMeshChunk',action='read')
    read(90,'(a)') line
    read(90,*) ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
    read(90,'(a)') line
    read(90,*) lon_center_chunk, lat_center_chunk, chunk_azi
    read(90,'(a)') line
    read(90,*) chunk_depth
    read(90,'(a)') line
    read(90,*) nel_lon,nel_lat, nel_depth
    read(90,'(a)') line
    read(90,'(a)') model1D_file
    read(90,'(a)') line
    read(90,*) buried_box
    if (buried_box) then
      read(90,'(a)') line
      read(90,*) radius_of_box_top
       radius_of_box_top =  radius_of_box_top * 1000.
    else
      radius_of_box_top = 6371000.
    endif
    model1D_file = 'MESH/'//trim(model1D_file)
    close(90)

    ! read 1D AxiSEM model
    call Read_dsm_model(model1D_file,vpv,vsv,density,zlayer,nlayer)

    ! modele 1D
    open(88,file='MESH/model_1D.in')
    write(88,*) nlayer,4
    do i=1,nlayer
      write(88,*) zlayer(i)
      write(88,'(4f20.10)') vpv(i,:)
      write(88,'(4f20.10)') vsv(i,:)
      write(88,'(4f20.10)') density(i,:)
    enddo
    z_bottom = minval(zgrid(:,:,:,:))
    write(88,*)  radius_of_box_top + z_bottom!6371000.+z_bottom
    write(88,*)  lon_center_chunk,  lat_center_chunk,  chunk_azi
    close(88)

    ! compute rotation matrix
    call compute_rotation_matrix(rotation_matrix,lon_center_chunk,lat_center_chunk, chunk_azi)

    open(91, file = 'MESH/list_ggl_boundary_spherical.txt')
    open(92, file = 'MESH/list_ggl_boundary_Cartesian.txt')
    open(89, file = 'MESH/flags_boundary.txt')

    open(90, file = 'MESH/Nb_ielm_faces.txt')
    write(90,*)  nspec2D_xmin
    write(90,*)  nspec2D_xmax
    write(90,*)  nspec2D_ymin
    write(90,*)  nspec2D_ymax
    write(90,*)  nspec2D_bottom
    write(90,*)  nspec2D_top
    close(90)

    ! xmin
    do ielm=1,nspec2D_xmin

      ispec=ibelm_xmin(ielm)

      write(89,*) ispec,ielm,1

      xelm(1)=xgrid(1,1,1,ispec)
      xelm(2)=xgrid(2,1,1,ispec)
      xelm(3)=xgrid(2,2,1,ispec)
      xelm(4)=xgrid(1,2,1,ispec)
      xelm(5)=xgrid(1,1,2,ispec)
      xelm(6)=xgrid(2,1,2,ispec)
      xelm(7)=xgrid(2,2,2,ispec)
      xelm(8)=xgrid(1,2,2,ispec)

      yelm(1)=ygrid(1,1,1,ispec)
      yelm(2)=ygrid(2,1,1,ispec)
      yelm(3)=ygrid(2,2,1,ispec)
      yelm(4)=ygrid(1,2,1,ispec)
      yelm(5)=ygrid(1,1,2,ispec)
      yelm(6)=ygrid(2,1,2,ispec)
      yelm(7)=ygrid(2,2,2,ispec)
      yelm(8)=ygrid(1,2,2,ispec)

      zelm(1)=zgrid(1,1,1,ispec)
      zelm(2)=zgrid(2,1,1,ispec)
      zelm(3)=zgrid(2,2,1,ispec)
      zelm(4)=zgrid(1,2,1,ispec)
      zelm(5)=zgrid(1,1,2,ispec)
      zelm(6)=zgrid(2,1,2,ispec)
      zelm(7)=zgrid(2,2,2,ispec)
      zelm(8)=zgrid(1,2,2,ispec)

      call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
      zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top !6371000.
      call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
      zstore(:,:,:) = zstore(:,:,:) -radius_of_box_top ! 6371000.
      call find_layer_in_axisem_model(ilayer,updown,radius(3,3,:),zlayer,nlayer)

      imin = 1
      imax = 1
      jmin = 1
      jmax = NGLLY
      kmin = 1
      kmax = NGLLZ

      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,1, &
                    ilayer,updown(k)
               write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
            enddo
         enddo
      enddo
    enddo

    ! xmax
    do ielm=1,nspec2D_xmax

      ispec=ibelm_xmax(ielm)

      write(89,*) ispec,ielm,2

      xelm(1)=xgrid(1,1,1,ispec)
      xelm(2)=xgrid(2,1,1,ispec)
      xelm(3)=xgrid(2,2,1,ispec)
      xelm(4)=xgrid(1,2,1,ispec)
      xelm(5)=xgrid(1,1,2,ispec)
      xelm(6)=xgrid(2,1,2,ispec)
      xelm(7)=xgrid(2,2,2,ispec)
      xelm(8)=xgrid(1,2,2,ispec)

      yelm(1)=ygrid(1,1,1,ispec)
      yelm(2)=ygrid(2,1,1,ispec)
      yelm(3)=ygrid(2,2,1,ispec)
      yelm(4)=ygrid(1,2,1,ispec)
      yelm(5)=ygrid(1,1,2,ispec)
      yelm(6)=ygrid(2,1,2,ispec)
      yelm(7)=ygrid(2,2,2,ispec)
      yelm(8)=ygrid(1,2,2,ispec)

      zelm(1)=zgrid(1,1,1,ispec)
      zelm(2)=zgrid(2,1,1,ispec)
      zelm(3)=zgrid(2,2,1,ispec)
      zelm(4)=zgrid(1,2,1,ispec)
      zelm(5)=zgrid(1,1,2,ispec)
      zelm(6)=zgrid(2,1,2,ispec)
      zelm(7)=zgrid(2,2,2,ispec)
      zelm(8)=zgrid(1,2,2,ispec)

      call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
      zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top !6371000.
      call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
      zstore(:,:,:) = zstore(:,:,:) -radius_of_box_top ! 6371000.
      call find_layer_in_axisem_model(ilayer,updown,radius(3,3,:),zlayer,nlayer)

      imin = NGLLX
      imax = NGLLX
      jmin = 1
      jmax = NGLLY
      kmin = 1
      kmax = NGLLZ

      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,2, &
                    ilayer,updown(k)
               write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
            enddo
         enddo
      enddo
    enddo

    ! ymin
    do ielm=1,nspec2D_ymin

      ispec=ibelm_ymin(ielm)

      write(89,*) ispec,ielm,3

      xelm(1)=xgrid(1,1,1,ispec)
      xelm(2)=xgrid(2,1,1,ispec)
      xelm(3)=xgrid(2,2,1,ispec)
      xelm(4)=xgrid(1,2,1,ispec)
      xelm(5)=xgrid(1,1,2,ispec)
      xelm(6)=xgrid(2,1,2,ispec)
      xelm(7)=xgrid(2,2,2,ispec)
      xelm(8)=xgrid(1,2,2,ispec)

      yelm(1)=ygrid(1,1,1,ispec)
      yelm(2)=ygrid(2,1,1,ispec)
      yelm(3)=ygrid(2,2,1,ispec)
      yelm(4)=ygrid(1,2,1,ispec)
      yelm(5)=ygrid(1,1,2,ispec)
      yelm(6)=ygrid(2,1,2,ispec)
      yelm(7)=ygrid(2,2,2,ispec)
      yelm(8)=ygrid(1,2,2,ispec)

      zelm(1)=zgrid(1,1,1,ispec)
      zelm(2)=zgrid(2,1,1,ispec)
      zelm(3)=zgrid(2,2,1,ispec)
      zelm(4)=zgrid(1,2,1,ispec)
      zelm(5)=zgrid(1,1,2,ispec)
      zelm(6)=zgrid(2,1,2,ispec)
      zelm(7)=zgrid(2,2,2,ispec)
      zelm(8)=zgrid(1,2,2,ispec)

      call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
      zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top !6371000.
      call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
      zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top !6371000.
      call find_layer_in_axisem_model(ilayer,updown,radius(3,3,:),zlayer,nlayer)

      imin = 1
      imax = NGLLX
      jmin = 1
      jmax = 1
      kmin = 1
      kmax = NGLLZ

      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,3, &
                    ilayer,updown(k)
               write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
            enddo
         enddo
      enddo
    enddo

    ! ymax
    do ielm=1,nspec2D_ymax

      ispec=ibelm_ymax(ielm)

      write(89,*) ispec,ielm,4

      xelm(1)=xgrid(1,1,1,ispec)
      xelm(2)=xgrid(2,1,1,ispec)
      xelm(3)=xgrid(2,2,1,ispec)
      xelm(4)=xgrid(1,2,1,ispec)
      xelm(5)=xgrid(1,1,2,ispec)
      xelm(6)=xgrid(2,1,2,ispec)
      xelm(7)=xgrid(2,2,2,ispec)
      xelm(8)=xgrid(1,2,2,ispec)

      yelm(1)=ygrid(1,1,1,ispec)
      yelm(2)=ygrid(2,1,1,ispec)
      yelm(3)=ygrid(2,2,1,ispec)
      yelm(4)=ygrid(1,2,1,ispec)
      yelm(5)=ygrid(1,1,2,ispec)
      yelm(6)=ygrid(2,1,2,ispec)
      yelm(7)=ygrid(2,2,2,ispec)
      yelm(8)=ygrid(1,2,2,ispec)

      zelm(1)=zgrid(1,1,1,ispec)
      zelm(2)=zgrid(2,1,1,ispec)
      zelm(3)=zgrid(2,2,1,ispec)
      zelm(4)=zgrid(1,2,1,ispec)
      zelm(5)=zgrid(1,1,2,ispec)
      zelm(6)=zgrid(2,1,2,ispec)
      zelm(7)=zgrid(2,2,2,ispec)
      zelm(8)=zgrid(1,2,2,ispec)

      call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
      zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top !6371000.
      call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
      zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top !6371000.
      call find_layer_in_axisem_model(ilayer,updown,radius(3,3,:),zlayer,nlayer)

      imin = 1
      imax = NGLLX
      jmin = NGLLY
      jmax = NGLLY
      kmin = 1
      kmax = NGLLZ

      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,4, &
                    ilayer,updown(k)
               write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
            enddo
         enddo
      enddo
    enddo

    ! bottom
    do ielm=1,nspec2D_BOTTOM

      ispec=ibelm_bottom(ielm)

      write(89,*) ispec,ielm,5

      xelm(1)=xgrid(1,1,1,ispec)
      xelm(2)=xgrid(2,1,1,ispec)
      xelm(3)=xgrid(2,2,1,ispec)
      xelm(4)=xgrid(1,2,1,ispec)
      xelm(5)=xgrid(1,1,2,ispec)
      xelm(6)=xgrid(2,1,2,ispec)
      xelm(7)=xgrid(2,2,2,ispec)
      xelm(8)=xgrid(1,2,2,ispec)

      yelm(1)=ygrid(1,1,1,ispec)
      yelm(2)=ygrid(2,1,1,ispec)
      yelm(3)=ygrid(2,2,1,ispec)
      yelm(4)=ygrid(1,2,1,ispec)
      yelm(5)=ygrid(1,1,2,ispec)
      yelm(6)=ygrid(2,1,2,ispec)
      yelm(7)=ygrid(2,2,2,ispec)
      yelm(8)=ygrid(1,2,2,ispec)

      zelm(1)=zgrid(1,1,1,ispec)
      zelm(2)=zgrid(2,1,1,ispec)
      zelm(3)=zgrid(2,2,1,ispec)
      zelm(4)=zgrid(1,2,1,ispec)
      zelm(5)=zgrid(1,1,2,ispec)
      zelm(6)=zgrid(2,1,2,ispec)
      zelm(7)=zgrid(2,2,2,ispec)
      zelm(8)=zgrid(1,2,2,ispec)

      call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
      zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top ! 6371000.
      call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
      zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top ! 6371000.
      call find_layer_in_axisem_model(ilayer,updown,radius(3,3,:),zlayer,nlayer)

      imin = 1
      imax = NGLLX
      jmin = 1
      jmax = NGLLY
      kmin = 1
      kmax = 1

      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,5, &
                    ilayer,updown(k)
               write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
            enddo
         enddo
      enddo
    enddo

    if (buried_box) then
      ! top
      do ielm=1,nspec2D_TOP

         ispec=ibelm_top(ielm)

         write(89,*) ispec,ielm,6

         xelm(1)=xgrid(1,1,1,ispec)
         xelm(2)=xgrid(2,1,1,ispec)
         xelm(3)=xgrid(2,2,1,ispec)
         xelm(4)=xgrid(1,2,1,ispec)
         xelm(5)=xgrid(1,1,2,ispec)
         xelm(6)=xgrid(2,1,2,ispec)
         xelm(7)=xgrid(2,2,2,ispec)
         xelm(8)=xgrid(1,2,2,ispec)

         yelm(1)=ygrid(1,1,1,ispec)
         yelm(2)=ygrid(2,1,1,ispec)
         yelm(3)=ygrid(2,2,1,ispec)
         yelm(4)=ygrid(1,2,1,ispec)
         yelm(5)=ygrid(1,1,2,ispec)
         yelm(6)=ygrid(2,1,2,ispec)
         yelm(7)=ygrid(2,2,2,ispec)
         yelm(8)=ygrid(1,2,2,ispec)

         zelm(1)=zgrid(1,1,1,ispec)
         zelm(2)=zgrid(2,1,1,ispec)
         zelm(3)=zgrid(2,2,1,ispec)
         zelm(4)=zgrid(1,2,1,ispec)
         zelm(5)=zgrid(1,1,2,ispec)
         zelm(6)=zgrid(2,1,2,ispec)
         zelm(7)=zgrid(2,2,2,ispec)
         zelm(8)=zgrid(1,2,2,ispec)

         call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
         zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top !6371000.
         call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
         zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top !6371000.
         call find_layer_in_axisem_model(ilayer,updown,radius(3,3,:),zlayer,nlayer)

         imin = 1
         imax = NGLLX
         jmin = 1
         jmax = NGLLY
         kmin = NGLLZ
         kmax = NGLLZ

         do k=kmin,kmax
            do j=jmin,jmax
               do i=imin,imax
                  write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,6, &
                       ilayer,updown(k)
                  write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
               enddo
            enddo
         enddo
      enddo
    endif

    close(89)
    close(91)
    close(92)

    deallocate(shape3D,dershape3D)
  endif

  end subroutine save_output_mesh_files_for_coupled_model


