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


module decompose_mesh_SCOTCH

  use part_decompose_mesh_SCOTCH

  implicit none

  include 'scotchf.h'

! number of partitions
  integer :: nparts ! e.g. 4 for partitioning for 4 CPUs or 4 processes

! mesh arrays
  integer(long) :: nspec
  integer, dimension(:,:), allocatable  :: elmnts
  integer, dimension(:,:), allocatable  :: mat
  integer, dimension(:), allocatable  :: part
  
  integer :: nnodes
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
  integer  ::  ninterfaces
  integer  :: my_ninterface

  integer(long)  :: nsize           ! Max number of elements that contain the same node.
  integer  :: nb_edges

  integer  :: ispec, inode
  integer  :: ngnod
  integer  :: max_neighbour         ! Real maximum number of neighbours per element
  integer(long)  :: sup_neighbour   ! Majoration of the maximum number of neighbours per element

  integer  :: ipart, nnodes_loc, nspec_loc
  integer  :: num_elmnt, num_node, num_mat

  ! boundaries
  integer  :: ispec2D
  integer  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, nspec2D_bottom, nspec2D_top
  integer, dimension(:), allocatable :: ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top
  integer, dimension(:,:), allocatable :: nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin
  integer, dimension(:,:), allocatable :: nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top

  ! moho surface (optional)
  integer :: nspec2D_moho
  integer, dimension(:), allocatable :: ibelm_moho
  integer, dimension(:,:), allocatable :: nodes_ibelm_moho

  character(len=256)  :: prname

  logical, dimension(:), allocatable :: mask_nodes_elmnts
  integer, dimension(:), allocatable :: used_nodes_elmnts

  double precision, dimension(SCOTCH_GRAPHDIM)  :: scotchgraph
  double precision, dimension(SCOTCH_STRATDIM)  :: scotchstrat
  character(len=256), parameter :: scotch_strategy='b{job=t,map=t,poli=S,sep=h{pass=30}}'
  integer  :: ier,idummy

  !pll
  double precision , dimension(:,:), allocatable :: mat_prop
  integer :: count_def_mat,count_undef_mat,imat
  character (len=30), dimension(:,:), allocatable :: undef_mat_prop

! default mesh file directory
  character(len=256) :: localpath_name
  character(len=256) :: outputpath_name

  integer :: aniso_flag,idomain_id
  double precision :: vp,vs,rho,qmu
! poroelastic parameters read in a new file
  double precision :: rhos,rhof,phi,tort,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,eta,mufr

  integer, parameter :: IIN_database = 15

  contains

  !----------------------------------------------------------------------------------------------
  ! reads in mesh files
  !----------------------------------------------------------------------------------------------
  subroutine read_mesh_files
    implicit none
    character(len=256)  :: line
    logical :: use_poroelastic_file

  ! sets number of nodes per element
    ngnod = esize

  ! reads node coordinates
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/nodes_coords_file',&
          status='old', form='formatted', iostat = ier)
    if( ier /= 0 ) then
      print*,'could not open file:',localpath_name(1:len_trim(localpath_name))//'/nodes_coords_file'
      stop 'error file open'
    endif
    read(98,*) nnodes
    allocate(nodes_coords(3,nnodes),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_coords'
    do inode = 1, nnodes
    ! format: #id_node #x_coordinate #y_coordinate #z_coordinate
      read(98,*) num_node, nodes_coords(1,num_node), nodes_coords(2,num_node), nodes_coords(3,num_node)
    !if(num_node /= inode)  stop "ERROR : Invalid nodes_coords file."
    end do
    close(98)
    print*, 'total number of nodes: '
    print*, '  nnodes = ', nnodes

  ! reads mesh elements indexing
  !(CUBIT calls this the connectivity, guess in the sense that it connects with the points index in
  ! the global coordinate file "nodes_coords_file"; it doesn't tell you which point is connected with others)
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/mesh_file', &
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening mesh_file'
    read(98,*) nspec
    allocate(elmnts(esize,nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array elmnts'
    do ispec = 1, nspec
      ! format: # element_id  #id_node1 ... #id_node8

      ! note: be aware that here we can have different node ordering for a cube element;
      !          the ordering from Cubit files might not be consistent for multiple volumes, or uneven, unstructured grids
      !
      !          guess here it assumes that spectral elements ordering is like first
      !          at the bottom of the element, anticlock-wise, i.e.
      !             point 1 = (0,0,0), point 2 = (0,1,0), point 3 = (1,1,0), point 4 = (1,0,0)
      !          then top (positive z-direction) of element
      !             point 5 = (0,0,1), point 6 = (0,1,1), point 7 = (1,1,1), point 8 = (1,0,1)

      !read(98,*) num_elmnt, elmnts(5,num_elmnt), elmnts(1,num_elmnt),elmnts(4,num_elmnt), elmnts(8,num_elmnt), &
      !      elmnts(6,num_elmnt), elmnts(2,num_elmnt), elmnts(3,num_elmnt), elmnts(7,num_elmnt)

      read(98,*) num_elmnt, elmnts(1,num_elmnt), elmnts(2,num_elmnt),elmnts(3,num_elmnt), elmnts(4,num_elmnt), &
            elmnts(5,num_elmnt), elmnts(6,num_elmnt), elmnts(7,num_elmnt), elmnts(8,num_elmnt)

      if((num_elmnt > nspec) .or. (num_elmnt < 1) )  stop "ERROR : Invalid mesh file."


      !outputs info for each element to see ordering
      !print*,'ispec: ',ispec
      !print*,'  ',num_elmnt, elmnts(5,num_elmnt), elmnts(1,num_elmnt),elmnts(4,num_elmnt), elmnts(8,num_elmnt), &
      !      elmnts(6,num_elmnt), elmnts(2,num_elmnt), elmnts(3,num_elmnt), elmnts(7,num_elmnt)
      !print*,'elem:',num_elmnt
      !do i=1,8
      !  print*,' i ',i,'val :',elmnts(i,num_elmnt),&
      !    nodes_coords(1,elmnts(i,num_elmnt)),nodes_coords(2,elmnts(i,num_elmnt)),nodes_coords(3,elmnts(i,num_elmnt))
      !enddo
      !print*

    end do
    close(98)
    print*, 'total number of spectral elements:'
    print*, '  nspec = ', nspec

  ! reads material associations
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/materials_file', &
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening materials_file'
    allocate(mat(2,nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array mat'
    mat(:,:) = 0
    do ispec = 1, nspec
      ! format: # id_element #flag
      ! note: be aware that elements may not be sorted in materials_file
      read(98,*) num_mat,mat(1,num_mat) !mat(1,ispec)!, mat(2,ispec)
      if((num_mat > nspec) .or. (num_mat < 1) ) stop "ERROR : Invalid mat file."
    end do
    close(98)

  ! TODO:
  ! must be changed, if  mat(1,i) < 0  1 == interface , 2 == tomography
    mat(2,:) = 1

  ! reads material definitions
  !
  ! note: format of nummaterial_velocity_file must be
  !
  ! #(1)material_domain_id #(2)material_id  #(3)rho  #(4)vp   #(5)vs   #(6)Q_mu  #(7)anisotropy_flag
  !
  ! where
  !     material_domain_id : 1=acoustic / 2=elastic / 3=poroelastic
  !     material_id               : number of material/volume
  !     rho                           : density
  !     vp                             : P-velocity
  !     vs                             : S-velocity
  !     Q_mu                      : 0=no attenuation
  !     anisotropy_flag        : 0=no anisotropy/ 1,2,.. check with implementation in aniso_model.f90
  ! Note that when poroelastic material, this file is a dummy except for material_domain_id & material_id,
  ! and that poroelastic materials are actually read from nummaterial_poroelastic_file, because CUBIT
  ! cannot support more than 10 attributes
    count_def_mat = 0
    count_undef_mat = 0
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/nummaterial_velocity_file',&
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening nummaterial_velocity_file'

    ! note: format #material_domain_id #material_id #...
    read(98,*,iostat=ier) idummy,num_mat
    print *,'materials:'
    ! counts materials (defined/undefined)
    do while (ier == 0)
       print*, '  num_mat = ',num_mat
       if(num_mat > 0 ) then
          ! positive materials_id: velocity values will be defined
          count_def_mat = count_def_mat + 1
       else
          ! negative materials_id: undefined material properties yet
          count_undef_mat = count_undef_mat + 1
       end if
       read(98,*,iostat=ier) idummy,num_mat
    end do
    close(98)
    print*, '  defined = ',count_def_mat, 'undefined = ',count_undef_mat
    ! check with material flags
    if( count_def_mat > 0 .and. maxval(mat(1,:)) > count_def_mat ) then
      print*,'error material definitions:'
      print*,'  materials associated in materials_file:',maxval(mat(1,:))
      print*,'  bigger than defined materials in nummaterial_velocity_file:',count_def_mat
      stop 'error materials'
    endif
    allocate(mat_prop(16,count_def_mat),stat=ier)
    if( ier /= 0 ) stop 'error allocating array mat_prop'
    allocate(undef_mat_prop(6,count_undef_mat),stat=ier)
    if( ier /= 0 ) stop 'error allocating array undef_mat_prop'
    mat_prop(:,:) = 0.d0
    undef_mat_prop(:,:) = ''

    ! reads in defined material properties
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/nummaterial_velocity_file', &
          status='old', form='formatted', iostat=ier)
    if( ier /= 0 ) stop 'error opening nummaterial_velocity_file'

  ! modif to read poro parameters, added if loop on idomain_id
  ! note: format of nummaterial_poroelastic_file located in MESH must be
  !
  ! #(1)rhos,#(2)rhof,#(3)phi,#(4)tort,#(5)kxx,#(6)kxy,#(7)kxz,#(8)kyy,#(9)kyz,#(10)kzz,
  ! #(11)kappas,#(12)kappaf,#(13)kappafr,#(14)eta,#(15)mufr
  !
  ! where
  !     rhos, rhof : solid & fluid density
  !     phi : porosity
  !     tort : tortuosity
  !     k : permeability tensor
  !     kappas, kappaf, kappafr : solid, fluid and frame bulk moduli
  !     eta : fluid viscosity
  !     mufr : frame shear modulus
    open(unit=97, file=localpath_name(1:len_trim(localpath_name))//'/nummaterial_poroelastic_file', &
          status='old', form='formatted', iostat=ier)
    ! checks if we can use file
    if( ier /= 0 ) then
      use_poroelastic_file = .false.
      !stop 'error opening nummaterial_poroelastic_file'
    else
      use_poroelastic_file = .true.
      print*, '  poroelastic material file found'
    endif
    ier = 0

    ! note: entries in nummaterial_velocity_file can be an unsorted list of all
    !          defined materials (material_id > 0) and undefined materials (material_id < 0 )
    do imat=1,count_def_mat
       ! material definitions
       !
       ! format: note that we save the arguments in a slightly different order in mat_prop(:,:)
       !              #(6) material_domain_id #(0) material_id  #(1) rho #(2) vp #(3) vs #(4) Q_mu #(5) anisotropy_flag
       !
       !read(98,*) idomain_id,num_mat,rho,vp,vs,qmu,aniso_flag
       ! reads lines unti it reaches a defined material
       num_mat = -1
       do while( num_mat < 0 .and. ier == 0)
         read(98,'(A256)',iostat=ier) line
         read(line,*) idomain_id,num_mat
       enddo
       if( ier /= 0 ) stop 'error reading in defined materials in nummaterial_velocity_file'

       ! reads in defined material properties
       read(line,*) idomain_id,num_mat,rho,vp,vs,qmu,aniso_flag

       ! checks material_id bounds
       if(num_mat < 1 .or. num_mat > count_def_mat)  stop "ERROR : Invalid nummaterial_velocity_file file."

       if(idomain_id == 1 .or. idomain_id == 2) then
         ! material is elastic or acoustic

         !read(98,*) num_mat, mat_prop(1,num_mat),mat_prop(2,num_mat),&
         !           mat_prop(3,num_mat),mat_prop(4,num_mat),mat_prop(5,num_mat)
         mat_prop(1,num_mat) = rho
         mat_prop(2,num_mat) = vp
         mat_prop(3,num_mat) = vs
         mat_prop(4,num_mat) = qmu
         mat_prop(5,num_mat) = aniso_flag
         mat_prop(6,num_mat) = idomain_id

       else
         ! material is poroelastic
         if( use_poroelastic_file .eqv. .false. ) stop 'error poroelastic material requires nummaterial_poroelastic_file'

         read(97,*) rhos,rhof,phi,tort,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,eta,mufr
         mat_prop(1,num_mat) = rhos
         mat_prop(2,num_mat) = rhof
         mat_prop(3,num_mat) = phi
         mat_prop(4,num_mat) = tort
         mat_prop(5,num_mat) = eta
         mat_prop(6,num_mat) = idomain_id
         mat_prop(7,num_mat) = kxx
         mat_prop(8,num_mat) = kxy
         mat_prop(9,num_mat) = kxz
         mat_prop(10,num_mat) = kyy
         mat_prop(11,num_mat) = kyz
         mat_prop(12,num_mat) = kzz
         mat_prop(13,num_mat) = kappas
         mat_prop(14,num_mat) = kappaf
         mat_prop(15,num_mat) = kappafr
         mat_prop(16,num_mat) = mufr

       endif !if(idomain_id == 1 .or. idomain_id == 2)

    end do

    ! reads in undefined material properties
    rewind(98,iostat=ier) ! back to the beginning of the file
    do imat=1,count_undef_mat
       !  undefined materials: have to be listed in decreasing order of material_id (start with -1, -2, etc...)
       !  format:
       !   - for interfaces
       !    #material_domain_id #material_id(<0) #type_name (="interface")
       !     #material_id_for_material_below #material_id_for_material_above
       !        example:     2  -1 interface 1 2
       !   - for tomography models
       !    #material_domain_id #material_id (<0) #type_name (="tomography") #block_name
       !        example:     2  -1 tomography elastic tomography_model.xyz 1
       ! reads lines unti it reaches a defined material
       num_mat = 1
       do while( num_mat >= 0 .and. ier == 0 )
         read(98,'(A256)',iostat=ier) line
         read(line,*) idomain_id,num_mat
       enddo
       if( ier /= 0 ) stop 'error reading in undefined materials in nummaterial_velocity_file'

       ! checks if interface or tomography definition
       read(line,*) undef_mat_prop(6,imat),undef_mat_prop(1,imat),undef_mat_prop(2,imat)
       if( trim(undef_mat_prop(2,imat)) == 'interface' ) then
         ! line will have 5 arguments, e.g.: 2  -1 interface 1 2
         read(line,*) undef_mat_prop(6,imat),undef_mat_prop(1,imat),undef_mat_prop(2,imat),&
                     undef_mat_prop(3,imat),undef_mat_prop(4,imat)
         undef_mat_prop(5,imat) = "0" ! dummy value
       else if( trim(undef_mat_prop(2,imat)) == 'tomography' ) then
         ! line will have 6 arguments, e.g.: 2  -1 tomography elastic tomography_model.xyz 1
         read(line,*) undef_mat_prop(6,imat),undef_mat_prop(1,imat),undef_mat_prop(2,imat),&
                        undef_mat_prop(3,imat),undef_mat_prop(4,imat),undef_mat_prop(5,imat)
       else
         stop "ERROR: invalid line in nummaterial_velocity_file for undefined material"
       endif

       !read(98,*) undef_mat_prop(6,imat),undef_mat_prop(1,imat),undef_mat_prop(2,imat),&
       !                 undef_mat_prop(3,imat),undef_mat_prop(4,imat),undef_mat_prop(5,imat)

       ! debug output
       !print*,'properties:'
       !print*,undef_mat_prop(:,imat)
       !print*

       ! checks material_id
       read(undef_mat_prop(1,imat),*) num_mat
       !print *,'material_id: ',num_mat
       if(num_mat > 0 .or. -num_mat > count_undef_mat)  &
            stop "ERROR : Invalid nummaterial_velocity_file for undefined materials."
       if(num_mat /= -imat)  &
            stop "ERROR : Invalid material_id in nummaterial_velocity_file for undefined materials."

       ! checks interface: flag_down/flag_up
       if( trim(undef_mat_prop(2,imat)) == 'interface' ) then
         ! flag_down
         read( undef_mat_prop(3,imat),*) num_mat
         if( num_mat > 0 ) then
          ! must point to a defined material
          if( num_mat > count_def_mat) &
               stop "ERROR: invalid flag_down in interface definition in nummaterial_velocity_file"
         else
          ! must point to an undefined material
          if( -num_mat > count_undef_mat) &
               stop "ERROR: invalid flag_down in interface definition in nummaterial_velocity_file"
         endif
         ! flag_up
         read( undef_mat_prop(4,imat),*) num_mat
         if( num_mat > 0 ) then
          ! must point to a defined material
          if( num_mat > count_def_mat) &
               stop "ERROR: invalid flag_up in interface definition in nummaterial_velocity_file"
         else
          ! must point to an undefined material
          if( -num_mat > count_undef_mat) &
               stop "ERROR: invalid flag_up in interface definition in nummaterial_velocity_file"
         endif
       endif
    end do
    if( use_poroelastic_file ) close(97)
    close(98)


    ! TODO:
    ! must be changed, if  mat(1,i) < 0  1 == interface , 2 == tomography
    do ispec=1,nspec
      ! get material_id
      num_mat = mat(1,ispec)
      if( num_mat < 0 ) then
        ! finds undefined material property
        do imat=1,count_undef_mat
          if( -imat == num_mat ) then
            ! interface
            if( trim(undef_mat_prop(2,imat)) == 'interface' ) then
              mat(2,ispec) = 1
            ! tomography
            elseif( trim(undef_mat_prop(2,imat)) == 'tomography' ) then
              mat(2,ispec) = 2
            else
              ! shouldn't encounter this case
              stop "error undefined material: type name not recognized"
            endif
          endif
        enddo
      else
        ! ispec belongs to a defined material
        mat(2,ispec) = 0
      endif
    enddo

  ! reads in absorbing boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_xmin', &
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) then
      nspec2D_xmin = 0
    else
      read(98,*) nspec2D_xmin
    endif
    allocate(ibelm_xmin(nspec2D_xmin),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_xmin'
    allocate(nodes_ibelm_xmin(4,nspec2D_xmin),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_xmin'
    do ispec2D = 1,nspec2D_xmin
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      ! note: ordering for CUBIT seems such that the normal of the face points outward of the element the face belongs to;
      !         in other words, nodes are in increasing order such that when looking from within the element outwards,
      !         they are ordered clockwise
      !
      !          doesn't necessarily have to start on top-rear, then bottom-rear, bottom-front, and finally top-front i.e.:
      !          point 1 = (0,1,1), point 2 = (0,1,0), point 3 = (0,0,0), point 4 = (0,0,1)
      read(98,*) ibelm_xmin(ispec2D), nodes_ibelm_xmin(1,ispec2D), nodes_ibelm_xmin(2,ispec2D), &
            nodes_ibelm_xmin(3,ispec2D), nodes_ibelm_xmin(4,ispec2D)

      !outputs info for each element for check of ordering
      !print*,'ispec2d:',ispec2d
      !print*,'  xmin:', ibelm_xmin(ispec2D), nodes_ibelm_xmin(1,ispec2D), nodes_ibelm_xmin(2,ispec2D), &
      !      nodes_ibelm_xmin(3,ispec2D), nodes_ibelm_xmin(4,ispec2D)
      !do i=1,4
      !  print*,'i',i,'val:',ibelm_xmin(ispec2d),nodes_coords(1,nodes_ibelm_xmin(i,ispec2D)), &
      !      nodes_coords(2,nodes_ibelm_xmin(i,ispec2D)),nodes_coords(3,nodes_ibelm_xmin(i,ispec2D))
      !enddo
      !print*
    end do
    close(98)
    print*, 'absorbing boundaries:'
    print*, '  nspec2D_xmin = ', nspec2D_xmin

  ! reads in absorbing boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_xmax', &
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) then
      nspec2D_xmax = 0
    else
      read(98,*) nspec2D_xmax
    endif
    allocate(ibelm_xmax(nspec2D_xmax),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_xmax'
    allocate(nodes_ibelm_xmax(4,nspec2D_xmax),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_xmax'
    do ispec2D = 1,nspec2D_xmax
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_xmax(ispec2D), nodes_ibelm_xmax(1,ispec2D), nodes_ibelm_xmax(2,ispec2D), &
            nodes_ibelm_xmax(3,ispec2D), nodes_ibelm_xmax(4,ispec2D)
    end do
    close(98)
    print*, '  nspec2D_xmax = ', nspec2D_xmax

  ! reads in absorbing boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_ymin', &
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) then
      nspec2D_ymin = 0
    else
      read(98,*) nspec2D_ymin
    endif
    allocate(ibelm_ymin(nspec2D_ymin),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_ymin'
    allocate(nodes_ibelm_ymin(4,nspec2D_ymin),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_ymin'
    do ispec2D = 1,nspec2D_ymin
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_ymin(ispec2D), nodes_ibelm_ymin(1,ispec2D), nodes_ibelm_ymin(2,ispec2D),  &
            nodes_ibelm_ymin(3,ispec2D), nodes_ibelm_ymin(4,ispec2D)
    end do
    close(98)
    print*, '  nspec2D_ymin = ', nspec2D_ymin

  ! reads in absorbing boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_ymax', &
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) then
      nspec2D_ymax = 0
    else
      read(98,*) nspec2D_ymax
    endif
    allocate(ibelm_ymax(nspec2D_ymax),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_ymax'
    allocate(nodes_ibelm_ymax(4,nspec2D_ymax),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_ymax'
    do ispec2D = 1,nspec2D_ymax
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_ymax(ispec2D), nodes_ibelm_ymax(1,ispec2D), nodes_ibelm_ymax(2,ispec2D),  &
            nodes_ibelm_ymax(3,ispec2D), nodes_ibelm_ymax(4,ispec2D)
    end do
    close(98)
    print*, '  nspec2D_ymax = ', nspec2D_ymax

  ! reads in absorbing boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_bottom', &
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) then
      nspec2D_bottom = 0
    else
      read(98,*) nspec2D_bottom
    endif
    allocate(ibelm_bottom(nspec2D_bottom),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_bottom'
    allocate(nodes_ibelm_bottom(4,nspec2D_bottom),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_bottom'
    do ispec2D = 1,nspec2D_bottom
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_bottom(ispec2D), nodes_ibelm_bottom(1,ispec2D), nodes_ibelm_bottom(2,ispec2D), &
            nodes_ibelm_bottom(3,ispec2D), nodes_ibelm_bottom(4,ispec2D)
    end do
    close(98)
    print*, '  nspec2D_bottom = ', nspec2D_bottom

  ! reads in free_surface boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/free_surface_file', &
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) then
      nspec2D_top = 0
    else
      read(98,*) nspec2D_top
    endif
    allocate(ibelm_top(nspec2D_top),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_top'
    allocate(nodes_ibelm_top(4,nspec2D_top),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_top'
    do ispec2D = 1,nspec2D_top
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_top(ispec2D), nodes_ibelm_top(1,ispec2D), nodes_ibelm_top(2,ispec2D), &
             nodes_ibelm_top(3,ispec2D), nodes_ibelm_top(4,ispec2D)
    end do
    close(98)
    print*, '  nspec2D_top = ', nspec2D_top

  ! reads in moho_surface boundary files (optional)
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/moho_surface_file', &
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) then
      nspec2D_moho = 0
    else
      read(98,*) nspec2D_moho
    endif
    allocate(ibelm_moho(nspec2D_moho),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_moho'
    allocate(nodes_ibelm_moho(4,nspec2D_moho),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_moho'
    do ispec2D = 1,nspec2D_moho
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_moho(ispec2D), nodes_ibelm_moho(1,ispec2D), nodes_ibelm_moho(2,ispec2D), &
             nodes_ibelm_moho(3,ispec2D), nodes_ibelm_moho(4,ispec2D)
    end do
    close(98)
    if( nspec2D_moho > 0 ) print*, '  nspec2D_moho = ', nspec2D_moho

  end subroutine read_mesh_files

  !----------------------------------------------------------------------------------------------
  ! checks valence of nodes
  !----------------------------------------------------------------------------------------------

  subroutine check_valence

    allocate(mask_nodes_elmnts(nnodes),stat=ier)
    if( ier /= 0 ) stop 'error allocating array mask_nodes_elmnts'
    allocate(used_nodes_elmnts(nnodes),stat=ier)
    if( ier /= 0 ) stop 'error allocating array used_nodes_elmnts'
    mask_nodes_elmnts(:) = .false.
    used_nodes_elmnts(:) = 0
    do ispec = 1, nspec
      do inode = 1, ESIZE
        mask_nodes_elmnts(elmnts(inode,ispec)) = .true.
        used_nodes_elmnts(elmnts(inode,ispec)) = used_nodes_elmnts(elmnts(inode,ispec)) + 1
      enddo
    enddo
    print *, 'nodes valence: '
    print *, '  min = ',minval(used_nodes_elmnts(:)),'max = ', maxval(used_nodes_elmnts(:))
    do inode = 1, nnodes
      if (.not. mask_nodes_elmnts(inode)) then
        stop 'ERROR : nodes not used.'
      endif
    enddo
    nsize = maxval(used_nodes_elmnts(:))
    sup_neighbour = ngnod * nsize - (ngnod + (ngnod/2 - 1)*nfaces)
    print*, '  nsize = ',nsize, 'sup_neighbour = ', sup_neighbour

  end subroutine check_valence

  !----------------------------------------------------------------------------------------------
  ! divides model into partitions using scotch library functions
  !----------------------------------------------------------------------------------------------

  subroutine scotch_partitioning

    implicit none
    ! local parameters
    integer, dimension(:),allocatable  :: num_material
    integer :: ier

    elmnts(:,:) = elmnts(:,:) - 1

    ! determines maximum neighbors based on 1 common node
    allocate(xadj(1:nspec+1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array xadj'
    allocate(adjncy(1:sup_neighbour*nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array adjncy'
    allocate(nnodes_elmnts(1:nnodes),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nnodes_elmnts'
    allocate(nodes_elmnts(1:nsize*nnodes),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_elmnts'

    call mesh2dual_ncommonnodes(nspec, nnodes, nsize, sup_neighbour, elmnts, xadj, adjncy, nnodes_elmnts, &
         nodes_elmnts, max_neighbour, 1)

    print*, 'mesh2dual: '
    print*, '  max_neighbour = ',max_neighbour

    nb_edges = xadj(nspec+1)

    ! allocates & initializes partioning of elements
    allocate(part(1:nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array part'
    part(:) = -1

    ! initializes
    ! elements load array
    allocate(elmnts_load(1:nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array elmnts_load'

    ! uniform load
    elmnts_load(:) = 1

    ! gets materials id associations
    allocate(num_material(1:nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array num_material'
    num_material(:) = mat(1,:)

    ! in case of acoustic/elastic/poro simulation, weights elements accordingly
    call acoustic_elastic_poro_load(elmnts_load,nspec,count_def_mat,count_undef_mat, &
                                  num_material,mat_prop,undef_mat_prop)

    deallocate(num_material)

    ! SCOTCH partitioning

    ! we use default strategy for partitioning, thus omit specifing explicit strategy .

    ! workflow preferred by F. Pellegrini (SCOTCH):
    !!This comes from the fact that, in version 5.1.8, the name
    !!for the "recursive bisection" method has changed from "b"
    !!("bipartitioning") to "r" ("recursive").
    !!
    !!As a general rule, do not try to set up strategies by
    !!yourself. The default strategy in Scotch will most probably
    !!provide better results. To use it, just call:
    !!
    !!SCOTCHFstratInit (),
    !!
    !!and use this "empty" strategy in the mapping routine
    !!(consequently, no call to SCOTCHFstratGraphMap () is
    !!required).
    !!
    !!This will make you independent from further changes (most
    !!probably improvements !;-)   ) in the strategy syntax.
    !!And you should see an improvement in performance, too,
    !!as your hand-made strategy did not make use of the
    !!multi-level framework.

    call scotchfstratinit (scotchstrat(1), ier)
     if (ier /= 0) then
       stop 'ERROR : MAIN : Cannot initialize strat'
    endif

    !call scotchfstratgraphmap (scotchstrat(1), trim(scotch_strategy), ier)
    ! if (ier /= 0) then
    !   stop 'ERROR : MAIN : Cannot build strat'
    !endif

    call scotchfgraphinit (scotchgraph (1), ier)
    if (ier /= 0) then
       stop 'ERROR : MAIN : Cannot initialize graph'
    endif

    ! fills graph structure : see user manual (scotch_user5.1.pdf, page 72/73)
    ! arguments: #(1) graph_structure       #(2) baseval(either 0/1)    #(3) number_of_vertices
    !                    #(4) adjacency_index_array         #(5) adjacency_end_index_array (optional)
    !                    #(6) vertex_load_array (optional) #(7) vertex_label_array
    !                    #(7) number_of_arcs                    #(8) adjacency_array
    !                    #(9) arc_load_array (optional)      #(10) ierror
    call scotchfgraphbuild (scotchgraph (1), 0, nspec, &
                          xadj (1), xadj (1), &
                          elmnts_load (1), xadj (1), &
                          nb_edges, adjncy (1), &
                          adjncy (1), ier)

    ! w/out element load, but adjacency array
    !call scotchfgraphbuild (scotchgraph (1), 0, nspec, &
    !                      xadj (1), xadj (1), &
    !                      xadj (1), xadj (1), &
    !                      nb_edges, adjncy (1), &
    !                      adjncy (1), ier)


    if (ier /= 0) then
       stop 'ERROR : MAIN : Cannot build graph'
    endif

    call scotchfgraphcheck (scotchgraph (1), ier)
    if (ier /= 0) then
       stop 'ERROR : MAIN : Invalid check'
    endif

    call scotchfgraphpart (scotchgraph (1), nparts, scotchstrat(1),part(1),ier)
    if (ier /= 0) then
       stop 'ERROR : MAIN : Cannot part graph'
    endif

    call scotchfgraphexit (scotchgraph (1), ier)
    if (ier /= 0) then
       stop 'ERROR : MAIN : Cannot destroy graph'
    endif

    call scotchfstratexit (scotchstrat(1), ier)
    if (ier /= 0) then
       stop 'ERROR : MAIN : Cannot destroy strat'
    endif
    
  ! re-partitioning puts poroelastic-elastic coupled elements into same partition
  !  integer  :: nfaces_coupled
  !  integer, dimension(:,:), pointer  :: faces_coupled
    call poro_elastic_repartitioning (nspec, nnodes, elmnts, &
                     count_def_mat, mat(1,:) , mat_prop, &
                     sup_neighbour, nsize, &
                     nparts, part)
                     !nparts, part, nfaces_coupled, faces_coupled)

  ! re-partitioning puts moho-surface coupled elements into same partition
    call moho_surface_repartitioning (nspec, nnodes, elmnts, &
                     sup_neighbour, nsize, nparts, part, &
                     nspec2D_moho,ibelm_moho,nodes_ibelm_moho )


  ! local number of each element for each partition    
    call build_glob2loc_elmnts(nspec, part, glob2loc_elmnts,nparts)

  ! local number of each node for each partition
    call build_glob2loc_nodes(nspec, nnodes,nsize, nnodes_elmnts, nodes_elmnts, part, &
         glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes, nparts)

  ! mpi interfaces
    ! acoustic/elastic/poroelastic boundaries will be split into different MPI partitions
    call build_interfaces(nspec, sup_neighbour, part, elmnts, &
                             xadj, adjncy, tab_interfaces, &
                             tab_size_interfaces, ninterfaces, &
                             nparts)

    !or: uncomment if you want acoustic/elastic boundaries NOT to be separated into different MPI partitions
    !call build_interfaces_no_ac_el_sep(nspec, sup_neighbour, part, elmnts, &
    !                          xadj, adjncy, tab_interfaces, &
    !                          tab_size_interfaces, ninterfaces, &
    !                          count_def_mat, mat_prop(3,:), mat(1,:), nparts)

  end subroutine scotch_partitioning

  !----------------------------------------------------------------------------------------------
  ! writes out new Databases files for each partition
  !----------------------------------------------------------------------------------------------

  subroutine write_mesh_databases

    implicit none
    !local parameters

    allocate(my_interfaces(0:ninterfaces-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array my_interfaces'
    allocate(my_nb_interfaces(0:ninterfaces-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array my_nb_interfaces'

    ! writes out Database file for each partition
    do ipart = 0, nparts-1

       ! opens output file
       write(prname, "(i6.6,'_Database')") ipart
       open(unit=IIN_database,file=outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname,&
            status='unknown', action='write', form='unformatted', iostat = ier)
       if( ier /= 0 ) then
        print*,'error file open:',outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname
        print*
        print*,'check if path exists:',outputpath_name(1:len_trim(outputpath_name))
        stop 'error file open Database'
       endif

       ! gets number of nodes
       call write_glob2loc_nodes_database(IIN_database, ipart, nnodes_loc, nodes_coords, &
                                  glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, nnodes, 1)

       ! gets number of spectral elements
       call write_partition_database(IIN_database, ipart, nspec_loc, nspec, elmnts, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, mat, ngnod, 1)

       ! writes out node coordinate locations
       !write(IIN_database,*) nnodes_loc
       write(IIN_database) nnodes_loc

       call write_glob2loc_nodes_database(IIN_database, ipart, nnodes_loc, nodes_coords,&
                                  glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, nnodes, 2)

       call write_material_props_database(IIN_database,count_def_mat,count_undef_mat, &
                                  mat_prop, undef_mat_prop)

       ! writes out spectral element indices
       !write(IIN_database,*) nspec_loc
       write(IIN_database) nspec_loc
       call write_partition_database(IIN_database, ipart, nspec_loc, nspec, elmnts, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, mat, ngnod, 2)

       ! writes out absorbing/free-surface boundaries
       call write_boundaries_database(IIN_database, ipart, nspec, nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, &
                                  nspec2D_ymax, nspec2D_bottom, nspec2D_top, &
                                  ibelm_xmin, ibelm_xmax, ibelm_ymin, &
                                  ibelm_ymax, ibelm_bottom, ibelm_top, &
                                  nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
                                  nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part)

       ! gets number of MPI interfaces
       call Write_interfaces_database(IIN_database, tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
                                  my_ninterface, my_interfaces, my_nb_interfaces, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, 1, nparts)

       ! writes out MPI interfaces elements
       !print*,' my interfaces:',my_ninterface,maxval(my_nb_interfaces)
       if( my_ninterface == 0 ) then
        !write(IIN_database,*) my_ninterface, 0
        write(IIN_database) my_ninterface, 0       ! avoids problem with maxval for empty array my_nb_interfaces
       else
        !write(IIN_database,*) my_ninterface, maxval(my_nb_interfaces)
        write(IIN_database) my_ninterface, maxval(my_nb_interfaces)
       endif

       call Write_interfaces_database(IIN_database, tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
                                  my_ninterface, my_interfaces, my_nb_interfaces, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, 2, nparts)

       ! writes out moho surface (optional)
       call write_moho_surface_database(IIN_database, ipart, nspec, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, &
                                  nspec2D_moho,ibelm_moho,nodes_ibelm_moho)

       close(IIN_database)

    end do
    print*, 'partitions: '
    print*, '  num = ',nparts
    print*
    print*, 'Databases files in directory: ',outputpath_name(1:len_trim(outputpath_name))
    print*, 'finished successfully'
    print*

  end subroutine write_mesh_databases

!end program pre_meshfem3D

end module decompose_mesh_SCOTCH

