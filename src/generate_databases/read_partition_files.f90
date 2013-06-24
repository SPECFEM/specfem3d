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
!


  subroutine read_partition_files

! reads in proc***_Databases files

  use generate_databases_par

  implicit none

  integer :: num_xmin,num_xmax,num_ymin,num_ymax,num_top,num_bottom,num
  integer :: num_cpml
  integer :: num_moho
  integer :: i,j
  !character(len=128) :: line

! read databases about external mesh simulation
! global node coordinates
  call create_name_database(prname,myrank,LOCAL_PATH)
  open(unit=IIN,file=prname(1:len_trim(prname))//'Database', &
        status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'rank ',myrank,' error opening file: ',prname(1:len_trim(prname))//'Database'
    print*,'please make sure file exists'
    call exit_mpi(myrank,'error opening database file')
  endif
  read(IIN) nnodes_ext_mesh

  allocate(nodes_coords_ext_mesh(NDIM,nnodes_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array nodes_coords_ext_mesh'

  do inode = 1, nnodes_ext_mesh
     read(IIN) dummy_node, nodes_coords_ext_mesh(1,inode), nodes_coords_ext_mesh(2,inode), &
                nodes_coords_ext_mesh(3,inode)
  enddo

  call sum_all_i(nnodes_ext_mesh,num)
  if(myrank == 0) then
    write(IMAIN,*) '  external mesh points: ',num
  endif
  call sync_all()

! read physical properties of the materials
! added poroelastic properties and filled with 0 the last 10 entries for elastic/acoustic
  read(IIN) nmat_ext_mesh, nundefMat_ext_mesh

  allocate(materials_ext_mesh(16,nmat_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array materials_ext_mesh'
  do imat = 1, nmat_ext_mesh
     ! (visco)elastic or acoustic format:
     ! #(1) rho   #(2) vp  #(3) vs  #(4) Q_mu  #(5) anisotropy_flag  #(6) material_domain_id  #(7) Q_kappa
     ! and remaining entries are filled with zeros.
     ! Q_kappa is not stored next to Q_mu for historical reasons, because it was added later.
     !
     ! poroelastic format:  rhos,rhof,phi,tort,eta,material_domain_id,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,mufr
     read(IIN) materials_ext_mesh(1,imat),  materials_ext_mesh(2,imat),  materials_ext_mesh(3,imat), &
          materials_ext_mesh(4,imat),  materials_ext_mesh(5,imat), materials_ext_mesh(6,imat), &
          materials_ext_mesh(7,imat),  materials_ext_mesh(8,imat), materials_ext_mesh(9,imat), &
          materials_ext_mesh(10,imat),  materials_ext_mesh(11,imat), materials_ext_mesh(12,imat), &
          materials_ext_mesh(13,imat),  materials_ext_mesh(14,imat), materials_ext_mesh(15,imat), &
          materials_ext_mesh(16,imat)
  enddo

  if(myrank == 0) then
    write(IMAIN,*) '  defined materials: ',nmat_ext_mesh
  endif
  call sync_all()

  allocate(undef_mat_prop(6,nundefMat_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array undef_mat_prop'
  do imat = 1, nundefMat_ext_mesh
     ! format example tomography:
     ! e.g.: -1 tomography elastic tomography_model.xyz 0 2
     ! format example interface:
     ! e.g.: -1 interface 14 15 0 2
     read(IIN) undef_mat_prop(1,imat),undef_mat_prop(2,imat),undef_mat_prop(3,imat),undef_mat_prop(4,imat), &
          undef_mat_prop(5,imat), undef_mat_prop(6,imat)
  enddo

  if(myrank == 0) then
    write(IMAIN,*) '  undefined materials: ',nundefMat_ext_mesh
  endif
  call sync_all()

! element indexing
  read(IIN) nelmnts_ext_mesh
  allocate(elmnts_ext_mesh(NGNOD,nelmnts_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array elmnts_ext_mesh'
  allocate(mat_ext_mesh(2,nelmnts_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array mat_ext_mesh'

  ! reads in material association for each spectral element and corner node indices
  do ispec = 1, nelmnts_ext_mesh
     ! format:
     ! # ispec_local # material_index_1 # material_index_2 # corner_id1 # corner_id2 # ... # corner_id8
     ! or
     ! # ispec_local # material_index_1 # material_index_2 # corner_id1 # corner_id2 # ... # corner_id27
     read(IIN) dummy_elmnt,mat_ext_mesh(1,ispec),mat_ext_mesh(2,ispec),(elmnts_ext_mesh(j,ispec),j=1,NGNOD)

     ! check debug
     if( dummy_elmnt /= ispec) stop "error ispec order in materials file"

  enddo
  NSPEC_AB = nelmnts_ext_mesh

  call sum_all_i(nspec_ab,num)
  if(myrank == 0) then
    write(IMAIN,*) ' total number of spectral elements: ',num
  endif
  call sync_all()

! reads absorbing/free-surface boundaries
  read(IIN) boundary_number ,nspec2D_xmin
  if(boundary_number /= 1) stop "Error : invalid database file"

  read(IIN) boundary_number ,nspec2D_xmax
  if(boundary_number /= 2) stop "Error : invalid database file"

  read(IIN) boundary_number ,nspec2D_ymin
  if(boundary_number /= 3) stop "Error : invalid database file"

  read(IIN) boundary_number ,nspec2D_ymax
  if(boundary_number /= 4) stop "Error : invalid database file"

  read(IIN) boundary_number ,nspec2D_bottom_ext
  if(boundary_number /= 5) stop "Error : invalid database file"

  read(IIN) boundary_number ,nspec2D_top_ext
  if(boundary_number /= 6) stop "Error : invalid database file"

  NSPEC2D_BOTTOM = nspec2D_bottom_ext
  NSPEC2D_TOP = nspec2D_top_ext

  allocate(ibelm_xmin(nspec2D_xmin),nodes_ibelm_xmin(NGNOD2D,nspec2D_xmin),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_xmin etc.'
  do ispec2D = 1,nspec2D_xmin
     read(IIN) ibelm_xmin(ispec2D),(nodes_ibelm_xmin(j,ispec2D),j=1,NGNOD2D)
  enddo

  allocate(ibelm_xmax(nspec2D_xmax),nodes_ibelm_xmax(NGNOD2D,nspec2D_xmax),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_xmax etc.'
  do ispec2D = 1,nspec2D_xmax
     read(IIN) ibelm_xmax(ispec2D),(nodes_ibelm_xmax(j,ispec2D),j=1,NGNOD2D)
  enddo

  allocate(ibelm_ymin(nspec2D_ymin),nodes_ibelm_ymin(NGNOD2D,nspec2D_ymin),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_ymin'
  do ispec2D = 1,nspec2D_ymin
     read(IIN) ibelm_ymin(ispec2D),(nodes_ibelm_ymin(j,ispec2D),j=1,NGNOD2D)
  enddo

  allocate(ibelm_ymax(nspec2D_ymax),nodes_ibelm_ymax(NGNOD2D,nspec2D_ymax),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_ymax etc.'
  do ispec2D = 1,nspec2D_ymax
     read(IIN) ibelm_ymax(ispec2D),(nodes_ibelm_ymax(j,ispec2D),j=1,NGNOD2D)
  enddo

  allocate(ibelm_bottom(nspec2D_bottom_ext),nodes_ibelm_bottom(NGNOD2D,nspec2D_bottom_ext),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_bottom etc.'
  do ispec2D = 1,nspec2D_bottom_ext
     read(IIN) ibelm_bottom(ispec2D),(nodes_ibelm_bottom(j,ispec2D),j=1,NGNOD2D)
  enddo

  allocate(ibelm_top(nspec2D_top_ext),nodes_ibelm_top(NGNOD2D,nspec2D_top_ext),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_top etc.'
  do ispec2D = 1,nspec2D_top_ext
     read(IIN) ibelm_top(ispec2D),(nodes_ibelm_top(j,ispec2D),j=1,NGNOD2D)
  enddo

  call sum_all_i(nspec2D_xmin,num_xmin)
  call sum_all_i(nspec2D_xmax,num_xmax)
  call sum_all_i(nspec2D_ymin,num_ymin)
  call sum_all_i(nspec2D_ymax,num_ymax)
  call sum_all_i(nspec2D_top_ext,num_top)
  call sum_all_i(nspec2D_bottom_ext,num_bottom)

  if(myrank == 0) then
    write(IMAIN,*) '  absorbing boundaries: '
    write(IMAIN,*) '    xmin,xmax: ',num_xmin,num_xmax
    write(IMAIN,*) '    ymin,ymax: ',num_ymin,num_ymax
    write(IMAIN,*) '    bottom,top: ',num_bottom,num_top
  endif
  call sync_all()

  ! reads number of C-PML elements in the global mesh
  nspec_cpml_tot = 0
  nspec_cpml = 0

  read(IIN) nspec_cpml_tot
  if(myrank == 0) then
     write(IMAIN,*) ' total number of C-PML elements in the global mesh: ',nspec_cpml_tot
  endif
  call sync_all()

  if( nspec_cpml_tot > 0 ) then
     ! reads number of C-PML elements in this partition
     read(IIN) nspec_cpml

     if(myrank == 0) then
        write(IMAIN,*) '  number of C-PML spectral elements in this partition: ',nspec_cpml
     endif
     call sync_all()

     call sum_all_i(nspec_cpml,num_cpml)

     ! checks that the sum of C-PML elements over all partitions is correct
     if( myrank == 0 .and. nspec_cpml_tot /= num_cpml ) stop 'error while summing C-PML elements over all partitions'

     ! reads C-PML regions and C-PML spectral elements global indexing
     allocate(CPML_to_spec(nspec_cpml),stat=ier)
     if(ier /= 0) stop 'error allocating array CPML_to_spec'
     allocate(CPML_regions(nspec_cpml),stat=ier)
     if(ier /= 0) stop 'error allocating array CPML_regions'

     do i=1,nspec_cpml
        ! #id_cpml_regions = 1 : X_surface C-PML
        ! #id_cpml_regions = 2 : Y_surface C-PML
        ! #id_cpml_regions = 3 : Z_surface C-PML
        ! #id_cpml_regions = 4 : XY_edge C-PML
        ! #id_cpml_regions = 5 : XZ_edge C-PML
        ! #id_cpml_regions = 6 : YZ_edge C-PML
        ! #id_cpml_regions = 7 : XYZ_corner C-PML
        !
        ! format: #id_cpml_element #id_cpml_regions
        read(IIN) CPML_to_spec(i), CPML_regions(i)
     enddo

     ! reads mask of C-PML elements for all elements in this partition
     allocate(is_CPML(NSPEC_AB),stat=ier)
     if(ier /= 0) stop 'error allocating array is_CPML'

     do i=1,NSPEC_AB
        read(IIN) is_CPML(i)
     enddo
  endif

! MPI interfaces between different partitions
  if( NPROC > 1 ) then
    ! format: #number_of_MPI_interfaces  #maximum_number_of_elements_on_each_interface
    read(IIN) num_interfaces_ext_mesh, max_interface_size_ext_mesh
  else
    num_interfaces_ext_mesh = 0
    max_interface_size_ext_mesh = 0
  endif

  ! allocates interfaces
  allocate(my_neighbours_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array my_neighbours_ext_mesh'
  allocate(my_nelmnts_neighbours_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array my_nelmnts_neighbours_ext_mesh'
  allocate(my_interfaces_ext_mesh(6,max_interface_size_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array my_interfaces_ext_mesh'
  allocate(ibool_interfaces_ext_mesh(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibool_interfaces_ext_mesh'
  allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array nibool_interfaces_ext_mesh'

  ! loops over MPI interfaces with other partitions
  do num_interface = 1, num_interfaces_ext_mesh
    ! format: #process_interface_id  #number_of_elements_on_interface
    ! where
    !     process_interface_id = rank of (neighbor) process to share MPI interface with
    !     number_of_elements_on_interface = number of interface elements
    read(IIN) my_neighbours_ext_mesh(num_interface), my_nelmnts_neighbours_ext_mesh(num_interface)

    ! loops over interface elements
    do ie = 1, my_nelmnts_neighbours_ext_mesh(num_interface)
      ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5)...
      !
      ! interface types:
      !     1  -  corner point only
      !     2  -  element edge
      !     4  -  element face
      read(IIN) my_interfaces_ext_mesh(1,ie,num_interface), my_interfaces_ext_mesh(2,ie,num_interface), &
                  my_interfaces_ext_mesh(3,ie,num_interface), my_interfaces_ext_mesh(4,ie,num_interface), &
                  my_interfaces_ext_mesh(5,ie,num_interface), my_interfaces_ext_mesh(6,ie,num_interface)
    enddo
  enddo

  call sum_all_i(num_interfaces_ext_mesh,num)
  if(myrank == 0) then
    write(IMAIN,*) '  number of MPI partition interfaces: ',num
  endif
  call sync_all()

  ! optional moho
  if( SAVE_MOHO_MESH ) then
    ! checks if additional line exists
    read(IIN,iostat=ier) boundary_number,nspec2D_moho_ext
    if( ier /= 0 ) then
      ! no moho informations given
      nspec2D_moho_ext = 0
      boundary_number = 7
    endif
    if(boundary_number /= 7) stop "Error : invalid database file"

    ! checks total number of elements
    call sum_all_i(nspec2D_moho_ext,num_moho)
    if( num_moho == 0 ) call exit_mpi(myrank,'error no moho mesh in database')

    ! reads in element informations
    allocate(ibelm_moho(nspec2D_moho_ext),nodes_ibelm_moho(NGNOD2D,nspec2D_moho_ext),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_moho etc.'
    do ispec2D = 1,nspec2D_moho_ext
      ! format: #element_id #node_id1 #node_id2 #node_id3 #node_id4
      read(IIN) ibelm_moho(ispec2D),(nodes_ibelm_moho(j,ispec2D),j=1,NGNOD2D)
    enddo

    ! user output
    if(myrank == 0) then
      write(IMAIN,*) '  moho surfaces: ',num_moho
    endif
    call sync_all()
  else
    ! allocate dummy array
    nspec2D_moho_ext = 0
    allocate(ibelm_moho(nspec2D_moho_ext),nodes_ibelm_moho(NGNOD2D,nspec2D_moho_ext),stat=ier)
    if( ier /= 0 ) stop 'error allocating dumy array ibelm_moho etc.'
  endif

  close(IIN)

  end subroutine read_partition_files
