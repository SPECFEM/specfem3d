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

  subroutine save_databases(prname,nspec,nglob,iproc_xi,iproc_eta, &
                            NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta,&
                            ibool,nodes_coords,ispec_material_id, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP,&
                            NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top,&
                            NMATERIALS,material_properties, &
                            nspec_CPML,CPML_to_spec,CPML_regions,is_CPML)

  use constants,only: MAX_STRING_LEN,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC

  implicit none

  include "constants_meshfem3D.h"

  ! number of spectral elements in each block
  integer nspec

  ! number of vertices in each block
  integer nglob

  ! MPI cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
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
  integer,dimension(2,nspec) :: material_index
  character(len=MAX_STRING_LEN), dimension(6,1) :: undef_mat_prop

  !------------------------------------------------------------------
  ! user parameter

  ! stores mesh files as cubit for single process run
  ! todo: we could put this parameter into the Mesh_Par_file
  logical, parameter :: SAVE_MESH_AS_CUBIT = .false.

  !------------------------------------------------------------------

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

  ! spectral-elements
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
     write(IIN_database) ibelm_xmin(i),ibool(1,1,1,ibelm_xmin(i)),ibool(1,NGLLY_M,1,ibelm_xmin(i)),&
          ibool(1,1,NGLLZ_M,ibelm_xmin(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_xmin(i))
  enddo
  do i=1,nspec2D_xmax
     write(IIN_database) ibelm_xmax(i),ibool(NGLLX_M,1,1,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,1,ibelm_xmax(i)), &
          ibool(NGLLX_M,1,NGLLZ_M,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_xmax(i))
  enddo
  do i=1,nspec2D_ymin
     write(IIN_database) ibelm_ymin(i),ibool(1,1,1,ibelm_ymin(i)),ibool(NGLLX_M,1,1,ibelm_ymin(i)),&
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
  if(nspec_CPML_total > 0) then
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
    ! determines number of mpi interfaces for each slice
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
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.))  then
          write(IIN_database) ispec,2,ibool(1,2,1,ispec),ibool(1,2,2,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(NE)) then
      write(IIN_database) addressing(iproc_xi+1,iproc_eta+1),nspec_interface(NE)
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.))  then
          write(IIN_database) ispec,2,ibool(2,2,1,ispec),ibool(2,2,2,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(SE)) then
      write(IIN_database) addressing(iproc_xi+1,iproc_eta-1),nspec_interface(SE)
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.))  then
          write(IIN_database) ispec,2,ibool(2,1,1,ispec),ibool(2,1,2,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(SW)) then
      write(IIN_database) addressing(iproc_xi-1,iproc_eta-1),nspec_interface(SW)
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.))  then
          write(IIN_database) ispec,2,ibool(1,1,1,ispec),ibool(1,1,2,ispec),-1,-1
        endif
      enddo
    endif

  else

    ! only 1 single slice, no mpi interfaces
    write(IIN_database) 0,0

    !! VM VM add outputs as CUBIT
    if (SAVE_MESH_AS_CUBIT) then
      call save_output_mesh_files_as_cubit(nspec,nglob, ibool,nodes_coords, ispec_material_id, &
                                           nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NMATERIALS,material_properties, &
                                           ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top)
    endif
  endif

  close(IIN_database)

  end subroutine save_databases

!---------------------------------------------------------------------------------------------------------------

  !! VM VM add subroutine for saving meshes in case of 1 mpi process
  subroutine save_output_mesh_files_as_cubit(nspec,nglob, ibool,nodes_coords, ispec_material_id, &
                                     nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                     NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NMATERIALS,material_properties, &
                                     ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top)

    use constants,only: MAX_STRING_LEN,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC

    implicit none

    include "constants_meshfem3D.h"

    integer, parameter :: IIN_database = 15
    ! number of spectral elements in each block
    integer :: nspec

    ! number of vertices in each block
    integer :: nglob

    ! MPI cartesian topology
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

    integer :: i,ispec,iglob,ier

    ! name of the database files

    integer :: mat_id,domain_id

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
         write(*,*) 'STOP : undefined mat not yet implemented '
         stop
       endif
    enddo
    close(IIN_database)

    open(IIN_database,file='MESH/materials_file')
    do ispec = 1, nspec
       write(IIN_database,*) ispec,ispec_material_id(ispec)
    enddo
    close(IIN_database)

    open(IIN_database,file='MESH/nodes_coords_file')
    write(IIN_database,*) nglob
    do iglob=1,nglob
       write(IIN_database,'(i14,3x,3(f20.5,1x))') iglob,nodes_coords(iglob,1),nodes_coords(iglob,2),nodes_coords(iglob,3)
    enddo
    close(IIN_database)

    open(IIN_database,file='MESH/mesh_file')
    write(IIN_database,*) nspec
    do ispec=1,nspec
       write(IIN_database,'(9i15)')  ispec,ibool(1,1,1,ispec),ibool(2,1,1,ispec),&
            ibool(2,2,1,ispec),ibool(1,2,1,ispec),&
            ibool(1,1,2,ispec),ibool(2,1,2,ispec),&
            ibool(2,2,2,ispec),ibool(1,2,2,ispec)
    enddo
    close(IIN_database)

    open(IIN_database,file='MESH/absorbing_surface_file_xmin')
    write(IIN_database,*) nspec2D_xmin
    do i=1,nspec2D_xmin
       write(IIN_database,'(5(i10,1x))') ibelm_xmin(i),ibool(1,1,1,ibelm_xmin(i)),ibool(1,NGLLY_M,1,ibelm_xmin(i)),&
          ibool(1,1,NGLLZ_M,ibelm_xmin(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_xmin(i))
    enddo
    close(IIN_database)

    open(IIN_database,file='MESH/absorbing_surface_file_xmax')
    write(IIN_database,*) nspec2D_xmax
    do i=1,nspec2D_xmax
       write(IIN_database,'(5(i10,1x))') ibelm_xmax(i),ibool(NGLLX_M,1,1,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,1,ibelm_xmax(i)), &
            ibool(NGLLX_M,1,NGLLZ_M,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_xmax(i))
    enddo
    close(IIN_database)

    open(IIN_database,file='MESH/absorbing_surface_file_ymin')
    write(IIN_database,*) nspec2D_ymin
    do i=1,nspec2D_ymin
       write(IIN_database,'(5(i10,1x))') ibelm_ymin(i),ibool(1,1,1,ibelm_ymin(i)),ibool(NGLLX_M,1,1,ibelm_ymin(i)),&
            ibool(1,1,NGLLZ_M,ibelm_ymin(i)),ibool(NGLLX_M,1,NGLLZ_M,ibelm_ymin(i))
    enddo
    close(IIN_database)

    open(IIN_database,file='MESH/absorbing_surface_file_ymax')
    write(IIN_database,*) nspec2D_ymax
    do i=1,nspec2D_ymax
       write(IIN_database,'(5(i10,1x))') ibelm_ymax(i),ibool(NGLLX_M,NGLLY_M,1,ibelm_ymax(i)),ibool(1,NGLLY_M,1,ibelm_ymax(i)), &
            ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_ymax(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_ymax(i))
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
       write(IIN_database,'(5(i10,1x))') ibelm_top(i),ibool(1,1,NGLLZ_M,ibelm_top(i)),ibool(NGLLX_M,1,NGLLZ_M,ibelm_top(i)), &
            ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_top(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_top(i))
    enddo
    close(IIN_database)

  end subroutine save_output_mesh_files_as_cubit
