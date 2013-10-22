!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

  subroutine save_arrays_solver_ext_mesh(nspec,nglob,APPROXIMATE_OCEAN_LOAD,ibool, &
                    num_interfaces_ext_mesh,my_neighbours_ext_mesh,nibool_interfaces_ext_mesh, &
                    max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                    SAVE_MESH_FILES,ANISOTROPY)

  use generate_databases_par, only: nspec_cpml,CPML_width_x,CPML_width_y,CPML_width_z,CPML_to_spec,&
                                    CPML_regions,is_CPML,nspec_cpml_tot, &
                                    d_store_x,d_store_y,d_store_z,k_store_x,k_store_y,k_store_z,&
                                    alpha_store_x,alpha_store_y,alpha_store_z, &
                                    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top,PML_CONDITIONS,&
                                    !for adjoint tomography
                                    SIMULATION_TYPE,SAVE_FORWARD,mask_ibool_interior_domain, &
                                    nglob_interface_PML_acoustic,points_interface_PML_acoustic,&
                                    nglob_interface_PML_elastic,points_interface_PML_elastic, &
                                    STACEY_ABSORBING_CONDITIONS
  use create_regions_mesh_ext_par

  implicit none

  integer :: nspec,nglob
  ! ocean load
  logical :: APPROXIMATE_OCEAN_LOAD
  ! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  ! MPI interfaces
  integer :: num_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer :: max_interface_size_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh

  logical :: SAVE_MESH_FILES
  logical :: ANISOTROPY

  ! local parameters
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer :: max_nibool_interfaces_ext_mesh

  integer :: ier,i
  character(len=256) :: filename

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
  if( ACOUSTIC_SIMULATION ) then
    write(IOUT) rmass_acoustic
  endif

! this array is needed for acoustic simulations but also for elastic simulations with CPML,
! thus we allocate it and read it in all cases (whether the simulation is acoustic, elastic, or acoustic/elastic)
  write(IOUT) rhostore

! elastic
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
  if( POROELASTIC_SIMULATION ) then
    write(IOUT) rmass_solid_poroelastic
    write(IOUT) rmass_fluid_poroelastic
    write(IOUT) rhoarraystore
    write(IOUT) kappaarraystore
    write(IOUT) etastore
    write(IOUT) tortstore
    write(IOUT) permstore
    write(IOUT) phistore
    write(IOUT) rho_vpI
    write(IOUT) rho_vpII
    write(IOUT) rho_vsI
  endif

! C-PML absorbing boundary conditions
  if( PML_CONDITIONS ) then
    write(IOUT) nspec_cpml
    write(IOUT) CPML_width_x
    write(IOUT) CPML_width_y
    write(IOUT) CPML_width_z
    if( nspec_cpml > 0 ) then
      write(IOUT) CPML_regions
      write(IOUT) CPML_to_spec
      write(IOUT) is_CPML
      write(IOUT) d_store_x
      write(IOUT) d_store_y
      write(IOUT) d_store_z
      write(IOUT) k_store_x
      write(IOUT) k_store_y
      write(IOUT) k_store_z
      write(IOUT) alpha_store_x
      write(IOUT) alpha_store_y
      write(IOUT) alpha_store_z
      ! --------------------------------------------------------------------------------------------
      ! for adjoint tomography
      ! save the array stored the points on interface between PML and interior computational domain
      ! --------------------------------------------------------------------------------------------
      if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        write(IOUT) nglob_interface_PML_acoustic
        write(IOUT) nglob_interface_PML_elastic
        if(nglob_interface_PML_acoustic > 0) write(IOUT) points_interface_PML_acoustic
        if(nglob_interface_PML_elastic > 0)  write(IOUT) points_interface_PML_elastic
      endif
    endif
  endif

! absorbing boundary surface
  write(IOUT) num_abs_boundary_faces
  if(PML_CONDITIONS)then
    if( num_abs_boundary_faces > 0 ) then
      write(IOUT) abs_boundary_ispec
      write(IOUT) abs_boundary_ijk
      write(IOUT) abs_boundary_jacobian2Dw
      write(IOUT) abs_boundary_normal
    endif
  else
    if( num_abs_boundary_faces > 0 ) then
      write(IOUT) abs_boundary_ispec
      write(IOUT) abs_boundary_ijk
      write(IOUT) abs_boundary_jacobian2Dw
      write(IOUT) abs_boundary_normal
      if( STACEY_ABSORBING_CONDITIONS ) then
        ! store mass matrix contributions
        if(ELASTIC_SIMULATION ) then
          write(IOUT) rmassx
          write(IOUT) rmassy
          write(IOUT) rmassz
        endif
        if(ACOUSTIC_SIMULATION) then
          write(IOUT) rmassz_acoustic
        endif
      endif
    endif
  endif

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
    write(IOUT) coupling_ac_el_ispec
    write(IOUT) coupling_ac_el_ijk
    write(IOUT) coupling_ac_el_jacobian2Dw
    write(IOUT) coupling_ac_el_normal
  endif

! acoustic-poroelastic coupling surface
  write(IOUT) num_coupling_ac_po_faces
  if( num_coupling_ac_po_faces > 0 ) then
    write(IOUT) coupling_ac_po_ispec
    write(IOUT) coupling_ac_po_ijk
    write(IOUT) coupling_ac_po_jacobian2Dw
    write(IOUT) coupling_ac_po_normal
  endif

! elastic-poroelastic coupling surface
  write(IOUT) num_coupling_el_po_faces
  if( num_coupling_el_po_faces > 0 ) then
    write(IOUT) coupling_el_po_ispec
    write(IOUT) coupling_po_el_ispec
    write(IOUT) coupling_el_po_ijk
    write(IOUT) coupling_po_el_ijk
    write(IOUT) coupling_el_po_jacobian2Dw
    write(IOUT) coupling_el_po_normal
  endif

  !MPI interfaces
  max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh(:))

  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array'
  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,i) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,i)
  enddo

  write(IOUT) num_interfaces_ext_mesh
  if( num_interfaces_ext_mesh > 0 ) then
    write(IOUT) max_nibool_interfaces_ext_mesh
    write(IOUT) my_neighbours_ext_mesh
    write(IOUT) nibool_interfaces_ext_mesh
    write(IOUT) ibool_interfaces_ext_mesh_dummy
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
    call save_arrays_solver_files(nspec,nglob,ibool)

    ! debug: saves 1. MPI interface
    !if( num_interfaces_ext_mesh >= 1 ) then
    !  filename = prname(1:len_trim(prname))//'MPI_1_points'
    !  call write_VTK_data_points(nglob, &
    !                    xstore_dummy,ystore_dummy,zstore_dummy, &
    !                    ibool_interfaces_ext_mesh_dummy(1:nibool_interfaces_ext_mesh(1),1), &
    !                    nibool_interfaces_ext_mesh(1), &
    !                    filename)
    !endif
  endif

  ! cleanup
  deallocate(ibool_interfaces_ext_mesh_dummy,stat=ier)
  if( ier /= 0 ) stop 'error deallocating array ibool_interfaces_ext_mesh_dummy'

  if( nspec_cpml_tot > 0 ) then
     deallocate(CPML_to_spec,stat=ier); if( ier /= 0 ) stop 'error deallocating array CPML_to_spec'
     deallocate(CPML_regions,stat=ier); if( ier /= 0 ) stop 'error deallocating array CPML_regions'
     deallocate(is_CPML,stat=ier); if( ier /= 0 ) stop 'error deallocating array is_CPML'
  endif

  if( PML_CONDITIONS ) then
     deallocate(d_store_x,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_x'
     deallocate(d_store_y,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_y'
     deallocate(d_store_z,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_z'
     deallocate(k_store_x,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_x'
     deallocate(k_store_y,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_y'
     deallocate(k_store_z,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_z'
     deallocate(alpha_store_x,stat=ier); if( ier /= 0 ) stop 'error deallocating array alpha_store_x'
     deallocate(alpha_store_y,stat=ier); if( ier /= 0 ) stop 'error deallocating array alpha_store_y'
     deallocate(alpha_store_z,stat=ier); if( ier /= 0 ) stop 'error deallocating array alpha_store_z'
     if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
       deallocate(mask_ibool_interior_domain,stat=ier)
       if(ier /= 0) stop 'error deallocating array mask_ibool_interior_domain'

       if(nglob_interface_PML_acoustic > 0) then
         deallocate(points_interface_PML_acoustic,stat=ier)
         if( ier /= 0 ) stop 'error deallocating array points_interface_PML_acoustic'
       endif

       if(nglob_interface_PML_elastic > 0) then
         deallocate(points_interface_PML_elastic,stat=ier)
         if( ier /= 0 ) stop 'error deallocating array points_interface_PML_elastic'
       endif
     endif
  endif

  end subroutine save_arrays_solver_ext_mesh

!
!-------------------------------------------------------------------------------------------------
!
  subroutine save_arrays_solver_files(nspec,nglob,ibool)

  use generate_databases_par, only: myrank
  use create_regions_mesh_ext_par

  implicit none

  integer :: nspec,nglob
  ! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: v_tmp
  integer,dimension(:),allocatable :: v_tmp_i
  integer :: ier,i
  integer, dimension(:), allocatable :: iglob_tmp
  integer :: j,inum
  character(len=256) :: filename

  logical,parameter :: DEBUG = .false.

  if( myrank == 0) then
    write(IMAIN,*) '     saving mesh files for AVS, OpenDX, Paraview'
    call flush_IMAIN()
  endif

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

  deallocate(v_tmp)

  ! VTK file output
  if( DEBUG ) then

    call sync_all()
    if( myrank == 0) then
      write(IMAIN,*) '     saving debugging mesh files'
      call flush_IMAIN()
    endif

    ! saves free surface points
    if( num_free_surface_faces > 0 ) then
      ! saves free surface interface points
      allocate( iglob_tmp(NGLLSQUARE*num_free_surface_faces),stat=ier)
      if( ier /= 0 ) stop 'error allocating array iglob_tmp'
      inum = 0
      iglob_tmp(:) = 0
      do i=1,num_free_surface_faces
        do j=1,NGLLSQUARE
          inum = inum+1
          iglob_tmp(inum) = ibool(free_surface_ijk(1,j,i), &
                                  free_surface_ijk(2,j,i), &
                                  free_surface_ijk(3,j,i), &
                                  free_surface_ispec(i) )
        enddo
      enddo
      filename = prname(1:len_trim(prname))//'free_surface'
      call write_VTK_data_points(nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy, &
                        iglob_tmp,NGLLSQUARE*num_free_surface_faces, &
                        filename)

      deallocate(iglob_tmp)
    endif

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

      deallocate(iglob_tmp,v_tmp_i)
    endif

    ! acoustic-poroelastic domains
    if( ACOUSTIC_SIMULATION .and. POROELASTIC_SIMULATION ) then
      ! saves points on acoustic-poroelastic coupling interface
      allocate( iglob_tmp(NGLLSQUARE*num_coupling_ac_po_faces),stat=ier)
      if( ier /= 0 ) stop 'error allocating array iglob_tmp'
      inum = 0
      iglob_tmp(:) = 0
      do i=1,num_coupling_ac_po_faces
        do j=1,NGLLSQUARE
          inum = inum+1
          iglob_tmp(inum) = ibool(coupling_ac_po_ijk(1,j,i), &
                                  coupling_ac_po_ijk(2,j,i), &
                                  coupling_ac_po_ijk(3,j,i), &
                                  coupling_ac_po_ispec(i) )
        enddo
      enddo
      filename = prname(1:len_trim(prname))//'coupling_acoustic_poroelastic'
      call write_VTK_data_points(nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy, &
                        iglob_tmp,NGLLSQUARE*num_coupling_ac_po_faces, &
                        filename)

      ! saves acoustic/poroelastic flag
      allocate(v_tmp_i(nspec),stat=ier)
      if( ier /= 0 ) stop 'error allocating array v_tmp_i'
      do i=1,nspec
        if( ispec_is_acoustic(i) ) then
          v_tmp_i(i) = 1
        else if( ispec_is_poroelastic(i) ) then
          v_tmp_i(i) = 2
        else
          v_tmp_i(i) = 0
        endif
      enddo
      filename = prname(1:len_trim(prname))//'acoustic_poroelastic_flag'
      call write_VTK_data_elem_i(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        v_tmp_i,filename)

      deallocate(v_tmp_i,iglob_tmp)
    endif !if( ACOUSTIC_SIMULATION .and. POROELASTIC_SIMULATION )

    ! elastic-poroelastic domains
    if( ELASTIC_SIMULATION .and. POROELASTIC_SIMULATION ) then
      ! saves points on elastic-poroelastic coupling interface
      allocate( iglob_tmp(NGLLSQUARE*num_coupling_el_po_faces),stat=ier)
      if( ier /= 0 ) stop 'error allocating array iglob_tmp'
      inum = 0
      iglob_tmp(:) = 0
      do i=1,num_coupling_el_po_faces
        do j=1,NGLLSQUARE
          inum = inum+1
          iglob_tmp(inum) = ibool(coupling_el_po_ijk(1,j,i), &
                                  coupling_el_po_ijk(2,j,i), &
                                  coupling_el_po_ijk(3,j,i), &
                                  coupling_el_po_ispec(i) )
        enddo
      enddo
      filename = prname(1:len_trim(prname))//'coupling_elastic_poroelastic'
      call write_VTK_data_points(nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy, &
                        iglob_tmp,NGLLSQUARE*num_coupling_el_po_faces, &
                        filename)

      ! saves elastic/poroelastic flag
      allocate(v_tmp_i(nspec),stat=ier)
      if( ier /= 0 ) stop 'error allocating array v_tmp_i'
      do i=1,nspec
        if( ispec_is_elastic(i) ) then
          v_tmp_i(i) = 1
        else if( ispec_is_poroelastic(i) ) then
          v_tmp_i(i) = 2
        else
          v_tmp_i(i) = 0
        endif
      enddo
      filename = prname(1:len_trim(prname))//'elastic_poroelastic_flag'
      call write_VTK_data_elem_i(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        v_tmp_i,filename)

      deallocate(v_tmp_i,iglob_tmp)
    endif !if( ACOUSTIC_SIMULATION .and. POROELASTIC_SIMULATION

  endif ! DEBUG

  end subroutine save_arrays_solver_files
