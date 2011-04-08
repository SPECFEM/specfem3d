!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

  subroutine save_arrays_solver_ext_mesh(nspec,nglob, &
                    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                    gammaxstore,gammaystore,gammazstore, &
                    jacobianstore, rho_vp,rho_vs,qmu_attenuation_store, &
                    rhostore,kappastore,mustore, &
                    rmass,rmass_acoustic,rmass_solid_poroelastic,rmass_fluid_poroelastic, &
                    OCEANS,rmass_ocean_load,NGLOB_OCEAN,&
                    ibool, &
                    xstore_dummy,ystore_dummy,zstore_dummy, &
                    abs_boundary_normal,abs_boundary_jacobian2Dw, &
                    abs_boundary_ijk,abs_boundary_ispec, &
                    num_abs_boundary_faces, &
                    free_surface_normal,free_surface_jacobian2Dw, &
                    free_surface_ijk,free_surface_ispec, &
                    num_free_surface_faces, &
                    coupling_ac_el_normal,coupling_ac_el_jacobian2Dw, &
                    coupling_ac_el_ijk,coupling_ac_el_ispec, &
                    num_coupling_ac_el_faces, &
                    num_interfaces_ext_mesh,my_neighbours_ext_mesh,nibool_interfaces_ext_mesh, &
                    max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                    prname,SAVE_MESH_FILES, &
                    ANISOTROPY,NSPEC_ANISO, &
                    c11store,c12store,c13store,c14store,c15store,c16store, &
                    c22store,c23store,c24store,c25store,c26store,c33store, &
                    c34store,c35store,c36store,c44store,c45store,c46store, &
                    c55store,c56store,c66store, &
                    ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic)

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
  logical :: OCEANS
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

! acoustic-elastic coupling surface
  integer :: num_coupling_ac_el_faces
  real(kind=CUSTOM_REAL) :: coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces)
  real(kind=CUSTOM_REAL) :: coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces)
  integer :: coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces)
  integer :: coupling_ac_el_ispec(num_coupling_ac_el_faces)

! MPI interfaces
  integer :: num_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer :: max_interface_size_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh

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
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
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
    if( OCEANS) then
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
  endif

! absorbing boundary surface
  write(IOUT) num_abs_boundary_faces
  write(IOUT) abs_boundary_ispec
  write(IOUT) abs_boundary_ijk
  write(IOUT) abs_boundary_jacobian2Dw
  write(IOUT) abs_boundary_normal

! free surface
  write(IOUT) num_free_surface_faces
  write(IOUT) free_surface_ispec
  write(IOUT) free_surface_ijk
  write(IOUT) free_surface_jacobian2Dw
  write(IOUT) free_surface_normal

! acoustic-elastic coupling surface
  write(IOUT) num_coupling_ac_el_faces
  write(IOUT) coupling_ac_el_ispec
  write(IOUT) coupling_ac_el_ijk
  write(IOUT) coupling_ac_el_jacobian2Dw
  write(IOUT) coupling_ac_el_normal

!MPI interfaces
  write(IOUT) num_interfaces_ext_mesh
  write(IOUT) maxval(nibool_interfaces_ext_mesh(:))
  write(IOUT) my_neighbours_ext_mesh
  write(IOUT) nibool_interfaces_ext_mesh

  allocate(ibool_interfaces_ext_mesh_dummy(maxval(nibool_interfaces_ext_mesh(:)),num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array'

  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,i) = ibool_interfaces_ext_mesh(1:maxval(nibool_interfaces_ext_mesh(:)),i)
  enddo
  write(IOUT) ibool_interfaces_ext_mesh_dummy

! anisotropy
  if( ANISOTROPY ) then
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

  close(IOUT)


! stores arrays in binary files
  if( SAVE_MESH_FILES ) then

    ! mesh arrays used for example in combine_vol_data.f90
    !--- x coordinate
    open(unit=27,file=prname(1:len_trim(prname))//'x.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file x.bin'
    write(27) xstore_dummy
    close(27)

    !--- y coordinate
    open(unit=27,file=prname(1:len_trim(prname))//'y.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file y.bin'
    write(27) ystore_dummy
    close(27)

    !--- z coordinate
    open(unit=27,file=prname(1:len_trim(prname))//'z.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file z.bin'
    write(27) zstore_dummy
    close(27)

    ! ibool
    open(unit=27,file=prname(1:len_trim(prname))//'ibool.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file ibool.bin'
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
    open(unit=27,file=prname(1:len_trim(prname))//'vp.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file vp.bin'
    write(27) v_tmp
    close(27)

    ! VTK file output
    ! vp values
    filename = prname(1:len_trim(prname))//'vp'
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
    open(unit=27,file=prname(1:len_trim(prname))//'vs.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file vs.bin'
    write(27) v_tmp
    close(27)

    ! VTK file output
    ! vs values
    filename = prname(1:len_trim(prname))//'vs'
    call write_VTK_data_gll_cr(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        v_tmp,filename)

    ! outputs density model for check
    v_tmp = 0.0
    where( rho_vp /= 0._CUSTOM_REAL ) v_tmp = rho_vp**2 / (FOUR_THIRDS * mustore + kappastore)
    open(unit=27,file=prname(1:len_trim(prname))//'rho.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rho.bin'
    write(27) v_tmp
    close(27)

    ! VTK file output
    ! saves attenuation flag assigned on each gll point into a vtk file
    filename = prname(1:len_trim(prname))//'attenuation'
    call write_VTK_data_gll_cr(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        qmu_attenuation_store,filename)

    ! VTK file output
    ! acoustic-elastic domains
    if( ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION ) then
      ! saves points on acoustic-elastic coupling interface
      allocate( iglob_tmp(NGLLSQUARE*num_coupling_ac_el_faces),stat=ier)
      if( ier /= 0 ) stop 'error allocating array iglob_tmp'
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
      allocate(v_tmp_i(nspec),stat=ier)
      if( ier /= 0 ) stop 'error allocating array v_tmp_i'
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

    !! saves 1. MPI interface
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

! cleanup
  deallocate(ibool_interfaces_ext_mesh_dummy,stat=ier); if( ier /= 0 ) stop 'error deallocating array'


  end subroutine save_arrays_solver_ext_mesh
