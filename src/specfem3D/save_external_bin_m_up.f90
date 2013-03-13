!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
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


! for external mesh

! Federica Magnoni:
! save_external_bin_m_up (compared to save_arrays_solver_ext_mesh)
! reads max_nibool_interfaces_ext_mesh instead of max_interface_size_ext_mesh

  subroutine save_external_bin_m_up(nspec,nglob, &
                    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                    gammaxstore,gammaystore,gammazstore, &
                    jacobianstore, rho_vp,rho_vs,qmu_attenuation_store, &
                    rhostore,kappastore,mustore, &
                    rmass,rmass_acoustic,rmass_solid_poroelastic,rmass_fluid_poroelastic, &
                    APPROXIMATE_OCEAN_LOAD,rmass_ocean_load,NGLOB_OCEAN,&
                    ibool, &
                    xstore_dummy,ystore_dummy,zstore_dummy, &
                    abs_boundary_normal,abs_boundary_jacobian2Dw, &
                    abs_boundary_ijk,abs_boundary_ispec, &
                    num_abs_boundary_faces, &
                    free_surface_normal,free_surface_jacobian2Dw, &
                    free_surface_ijk,free_surface_ispec, &
                    num_free_surface_faces, &
                    num_interfaces_ext_mesh,my_neighbours_ext_mesh,nibool_interfaces_ext_mesh, &
                    max_nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                    prname,SAVE_MESH_FILES, &
                    ANISOTROPY,NSPEC_ANISO, &
                    c11store,c12store,c13store,c14store,c15store,c16store, &
                    c22store,c23store,c24store,c25store,c26store,c33store, &
                    c34store,c35store,c36store,c44store,c45store,c46store, &
                    c55store,c56store,c66store, &
                    ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic)

  use specfem_par,only: &
    ispec_is_inner

  use specfem_par_elastic,only: &
    rmassx,rmassy,rmassz, &
    nspec_inner_elastic,nspec_outer_elastic,num_phase_ispec_elastic,phase_ispec_inner_elastic, &
    num_colors_outer_elastic,num_colors_inner_elastic,num_elem_colors_elastic

  use specfem_par_acoustic,only: &
    rmassz_acoustic,num_coupling_ac_po_faces, &
    num_coupling_ac_el_faces,coupling_ac_el_ijk,coupling_ac_el_ispec, &
    nspec_inner_acoustic,nspec_outer_acoustic,num_phase_ispec_acoustic,phase_ispec_inner_acoustic, &
    num_colors_outer_acoustic,num_colors_inner_acoustic,num_elem_colors_acoustic

  use specfem_par_poroelastic,only: &
    num_coupling_el_po_faces, &
    nspec_inner_poroelastic,nspec_outer_poroelastic,num_phase_ispec_poroelastic,phase_ispec_inner_poroelastic

  implicit none

  include "constants.h"

  integer :: nspec,nglob

! jacobian
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xixstore,xiystore,xizstore, &
            etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rho_vp,rho_vs

! attenuation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: qmu_attenuation_store

! material
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rhostore,kappastore,mustore
  real(kind=CUSTOM_REAL), dimension(nglob) :: rmass,rmass_acoustic, &
            rmass_solid_poroelastic,rmass_fluid_poroelastic

! ocean load
  logical :: APPROXIMATE_OCEAN_LOAD
  integer :: NGLOB_OCEAN
  real(kind=CUSTOM_REAL),dimension(NGLOB_OCEAN) :: rmass_ocean_load

! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! absorbing boundary surface
  integer :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)

! free surface
  integer :: num_free_surface_faces
  real(kind=CUSTOM_REAL) :: free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces)
  real(kind=CUSTOM_REAL) :: free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ispec(num_free_surface_faces)

! MPI interfaces
  integer :: num_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer :: max_nibool_interfaces_ext_mesh  !magnoni
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh  !magnoni


! file name
  character(len=256) prname
  logical :: SAVE_MESH_FILES

! anisotropy
  logical :: ANISOTROPY
  integer :: NSPEC_ANISO
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO) :: &
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store, &
            c34store,c35store,c36store,c44store,c45store,c46store, &
            c55store,c56store,c66store

! material domain flags
  logical, dimension(nspec) :: ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic

! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: v_tmp
  integer,dimension(:),allocatable :: v_tmp_i

  !real(kind=CUSTOM_REAL) :: minimum(1)
  integer :: ier,i
  logical :: ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION
  character(len=256) :: filename

  integer, dimension(:), allocatable :: iglob_tmp
  integer :: j,inum

! saves mesh file proc***_external_mesh.bin
  filename = prname(1:len_trim(prname))//'external_mesh.bin'
  open(unit=IOUT,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
  if( ier /= 0 ) stop 'error opening database proc######_external_mesh.bin'

  write(IOUT) nspec
  write(IOUT) nglob

  write(IOUT) ibool

  write(IOUT) xstore_dummy
  write(IOUT) ystore_dummy
  write(IOUT) zstore_dummy

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

  write(IOUT) kappastore
  write(IOUT) mustore

  write(IOUT) ispec_is_acoustic
  write(IOUT) ispec_is_elastic
  write(IOUT) ispec_is_poroelastic

! acoustic
! all processes will have acoustic_simulation set if any flag is .true. somewhere
  call any_all_l( ANY(ispec_is_acoustic), ACOUSTIC_SIMULATION )
  if( ACOUSTIC_SIMULATION ) then
    write(IOUT) rmass_acoustic
    write(IOUT) rhostore
  endif

! elastic
  call any_all_l( ANY(ispec_is_elastic), ELASTIC_SIMULATION )
  if( ELASTIC_SIMULATION ) then
    write(IOUT) rmass
    if( APPROXIMATE_OCEAN_LOAD) then
      write(IOUT) rmass_ocean_load
    endif
    !pll Stacey
    write(IOUT) rho_vp
    write(IOUT) rho_vs
  endif

! poroelastic
  call any_all_l( ANY(ispec_is_poroelastic), POROELASTIC_SIMULATION )
  if( POROELASTIC_SIMULATION ) then
    write(IOUT) rmass_solid_poroelastic
    write(IOUT) rmass_fluid_poroelastic
    stop 'model update with poroelastic domains not supported yet'
  endif

! absorbing boundary surface
  write(IOUT) num_abs_boundary_faces
  if(num_abs_boundary_faces > 0 ) then
    write(IOUT) abs_boundary_ispec
    write(IOUT) abs_boundary_ijk
    write(IOUT) abs_boundary_jacobian2Dw
    write(IOUT) abs_boundary_normal
    ! store mass matrix contributions
    if(ELASTIC_SIMULATION) then
     write(IOUT) rmassx
     write(IOUT) rmassy
     write(IOUT) rmassz
    endif
    if(ACOUSTIC_SIMULATION) then
      write(IOUT) rmassz_acoustic
    endif
  endif

! free surface
  write(IOUT) num_free_surface_faces
  if( num_free_surface_faces > 0 ) then
    write(IOUT) free_surface_ispec
    write(IOUT) free_surface_ijk
    write(IOUT) free_surface_jacobian2Dw
    write(IOUT) free_surface_normal
  endif

! acoustic-elastic coupling surface
  write(IOUT) num_coupling_ac_el_faces
  if( num_coupling_ac_el_faces > 0 ) then
    stop 'coupling ac_po not updated yet'
  endif

! acoustic-poroelastic coupling surface
  write(IOUT) num_coupling_ac_po_faces
  if( num_coupling_ac_po_faces > 0 ) then
    stop 'coupling ac_po not updated yet'
  endif

! elastic-poroelastic coupling surface
  write(IOUT) num_coupling_el_po_faces
  if( num_coupling_el_po_faces > 0 ) then
    stop 'coupling ac_po not updated yet'
  endif

!MPI interfaces
  write(IOUT) num_interfaces_ext_mesh
  if( num_interfaces_ext_mesh > 0 ) then
    write(IOUT) max_nibool_interfaces_ext_mesh
    write(IOUT) my_neighbours_ext_mesh
    write(IOUT) nibool_interfaces_ext_mesh
    write(IOUT) ibool_interfaces_ext_mesh   !magnoni
  endif

! anisotropy
  if( ELASTIC_SIMULATION .and. ANISOTROPY ) then
    write(IOUT) c11store
    write(IOUT) c12store
    write(IOUT) c13store
    write(IOUT) c14store
    write(IOUT) c15store
    write(IOUT) c16store
    write(IOUT) c22store
    write(IOUT) c23store
    write(IOUT) c24store
    write(IOUT) c25store
    write(IOUT) c26store
    write(IOUT) c33store
    write(IOUT) c34store
    write(IOUT) c35store
    write(IOUT) c36store
    write(IOUT) c44store
    write(IOUT) c45store
    write(IOUT) c46store
    write(IOUT) c55store
    write(IOUT) c56store
    write(IOUT) c66store
  endif

! inner/outer elements
  write(IOUT) ispec_is_inner

  if( ACOUSTIC_SIMULATION ) then
    write(IOUT) nspec_inner_acoustic,nspec_outer_acoustic
    write(IOUT) num_phase_ispec_acoustic
    if(num_phase_ispec_acoustic > 0 ) write(IOUT) phase_ispec_inner_acoustic
  endif

  if( ELASTIC_SIMULATION ) then
    write(IOUT) nspec_inner_elastic,nspec_outer_elastic
    write(IOUT) num_phase_ispec_elastic
    if(num_phase_ispec_elastic > 0 ) write(IOUT) phase_ispec_inner_elastic
  endif

  if( POROELASTIC_SIMULATION ) then
    write(IOUT) nspec_inner_poroelastic,nspec_outer_poroelastic
    write(IOUT) num_phase_ispec_poroelastic
    if(num_phase_ispec_poroelastic > 0 ) write(IOUT) phase_ispec_inner_poroelastic
  endif

  ! mesh coloring
  if( USE_MESH_COLORING_GPU ) then
    if( ACOUSTIC_SIMULATION ) then
      write(IOUT) num_colors_outer_acoustic,num_colors_inner_acoustic
      write(IOUT) num_elem_colors_acoustic
    endif
    if( ELASTIC_SIMULATION ) then
      write(IOUT) num_colors_outer_elastic,num_colors_inner_elastic
      write(IOUT) num_elem_colors_elastic
    endif
  endif

  close(IOUT)


! stores arrays in binary files
  if( SAVE_MESH_FILES ) then

    ! mesh arrays used for example in combine_vol_data.f90
    !--- x coordinate
    open(unit=27,file=prname(1:len_trim(prname))//'x.bin',status='unknown',form='unformatted')
    write(27) xstore_dummy
    close(27)

    !--- y coordinate
    open(unit=27,file=prname(1:len_trim(prname))//'y.bin',status='unknown',form='unformatted')
    write(27) ystore_dummy
    close(27)

    !--- z coordinate
    open(unit=27,file=prname(1:len_trim(prname))//'z.bin',status='unknown',form='unformatted')
    write(27) zstore_dummy
    close(27)

    ! ibool
    open(unit=27,file=prname(1:len_trim(prname))//'ibool.bin',status='unknown',form='unformatted')
    write(27) ibool
    close(27)

    allocate( v_tmp(NGLLX,NGLLY,NGLLZ,nspec), stat=ier); if( ier /= 0 ) stop 'error allocating array '

    ! vp (for checking the mesh and model)
    !minimum = minval( abs(rho_vp) )
    !if( minimum(1) /= 0.0 ) then
    !  v_tmp = (FOUR_THIRDS * mustore + kappastore) / rho_vp
    !else
    !  v_tmp = 0.0
    !endif
    v_tmp = 0.0
    where( rho_vp /= 0._CUSTOM_REAL ) v_tmp = (FOUR_THIRDS * mustore + kappastore) / rho_vp
    open(unit=27,file=prname(1:len_trim(prname))//'vp_new.bin',status='unknown',form='unformatted')
    write(27) v_tmp
    close(27)

    ! VTK file output
    ! vp values
    filename = prname(1:len_trim(prname))//'vp_new'
    call write_VTK_data_gll_cr(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        v_tmp,filename)


    ! vs (for checking the mesh and model)
    !minimum = minval( abs(rho_vs) )
    !if( minimum(1) /= 0.0 ) then
    !  v_tmp = mustore / rho_vs
    !else
    !  v_tmp = 0.0
    !endif
    v_tmp = 0.0
    where( rho_vs /= 0._CUSTOM_REAL )  v_tmp = mustore / rho_vs
    open(unit=27,file=prname(1:len_trim(prname))//'vs_new.bin',status='unknown',form='unformatted')
    write(27) v_tmp
    close(27)

    ! VTK file output
    ! vs values
    filename = prname(1:len_trim(prname))//'vs_new'
    call write_VTK_data_gll_cr(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        v_tmp,filename)

    ! outputs density model for check
    v_tmp = 0.0
    where( rho_vp /= 0._CUSTOM_REAL ) v_tmp = rho_vp**2 / (FOUR_THIRDS * mustore + kappastore)
    open(unit=27,file=prname(1:len_trim(prname))//'rho_new.bin',status='unknown',form='unformatted')
    write(27) v_tmp
    close(27)

    ! VTK file output
    ! density model
    filename = prname(1:len_trim(prname))//'rho_new'
    call write_VTK_data_gll_cr(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        v_tmp,filename)


    ! VTK file output
    ! saves attenuation flag assigned on each gll point into a vtk file
    filename = prname(1:len_trim(prname))//'attenuation'
    call write_VTK_data_gll_cr(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        qmu_attenuation_store,filename)

!     !magnoni
!     !check
!
!       !jacobian
!       filename = prname(1:len_trim(prname))//'jacobian'
!       call write_VTK_data_gll_cr(nspec,nglob, &
!                           xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!                           jacobianstore,filename)
!
!     !check
!     !magnoni


    ! VTK file output
    ! acoustic-elastic domains
    if( ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION ) then
      ! saves points on acoustic-elastic coupling interface
      allocate( iglob_tmp(NGLLSQUARE*num_coupling_ac_el_faces))
      inum = 0
      iglob_tmp(:) = 0
      do i=1,num_coupling_ac_el_faces
        do j=1,NGLLSQUARE
          inum = inum+1
          iglob_tmp(inum) = ibool(coupling_ac_el_ijk(1,j,i), &
                                  coupling_ac_el_ijk(2,j,i), &
                                  coupling_ac_el_ijk(3,j,i), &
                                  coupling_ac_el_ispec(i) )
        enddo
      enddo
      filename = prname(1:len_trim(prname))//'coupling_acoustic_elastic'
      call write_VTK_data_points(nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy, &
                        iglob_tmp,NGLLSQUARE*num_coupling_ac_el_faces, &
                        filename)

      ! saves acoustic/elastic flag
      allocate(v_tmp_i(nspec))
      do i=1,nspec
        if( ispec_is_acoustic(i) ) then
          v_tmp_i(i) = 1
        else if( ispec_is_elastic(i) ) then
          v_tmp_i(i) = 2
        else
          v_tmp_i(i) = 0
        endif
      enddo
      filename = prname(1:len_trim(prname))//'acoustic_elastic_flag'
      call write_VTK_data_elem_i(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        v_tmp_i,filename)
    endif

    !debug: saves 1. MPI interface
    !    if( num_interfaces_ext_mesh >= 1 ) then
    !      filename = prname(1:len_trim(prname))//'MPI_1_points'
    !      call write_VTK_data_points(nglob, &
    !                        xstore_dummy,ystore_dummy,zstore_dummy, &
    !                        ibool_interfaces_ext_mesh_dummy(1:nibool_interfaces_ext_mesh(1),1), &
    !                        nibool_interfaces_ext_mesh(1), &
    !                        filename)
    !    endif
    !

    deallocate(v_tmp)

  endif ! SAVE_MESH_FILES


  end subroutine save_external_bin_m_up
