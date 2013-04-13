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

!! DK DK Oct 2012: this can significantly reduce the time it takes to decompose the mesh in the case of very large meshes,
!! DK DK Oct 2012: and in principle it is safe to do that; however it is not fully tested yet
!! DK DK Oct 2012: and thus I comment it out for now
!!!! #define USE_TWO_CALLS_TO_mesh2dual

!! DK DK added support for METIS as well, in addition to SCOTCH
!! DK DK thus uncomment the define statement below if you want to use METIS instead of SCOTCH
!! DK DK and then also add the path to the METIS library files in Makefile.in in this directory
!! DK DK before running the "configure" script
!!!! #define USE_METIS_INSTEAD_OF_SCOTCH

module decompose_mesh

  use part_decompose_mesh
  use fault_scotch

  implicit none

  include 'scotchf.h'

!! DK DK added this because poroelastic repartitioning routine of Christina Morency is currently broken
! implement mesh repartitioning of poroelastic-elastic interface
! (the risk being to break the nice load balancing created by the domain decomposer for high-performance computing)
  logical, parameter :: PORO_INTERFACE_REPARTITIONING = .false.

! number of partitions
  integer :: nparts

! mesh arrays
  integer :: nspec
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
  integer ::  ninterfaces
  integer :: my_ninterface

  integer :: nsize           ! max number of elements that contain the same node
  integer :: nb_edges

  integer :: ispec, inode
  integer :: max_neighbour   ! real maximum number of neighbours per element
  integer :: sup_neighbour   ! majoration (overestimate) of the maximum number of neighbours per element

  integer :: ipart, nnodes_loc, nspec_local,ncommonnodes
  integer :: num_elmnt, num_node, num_mat

  ! boundaries
  integer :: ispec2D
  integer :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, nspec2D_bottom, nspec2D_top
  integer, dimension(:), allocatable :: ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top
  integer, dimension(:,:), allocatable :: nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin
  integer, dimension(:,:), allocatable :: nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top

  ! C-PML absorbing boundary conditions
  integer :: ispec_CPML
  integer :: nspec_cpml
  integer, dimension(:), allocatable :: CPML_to_spec, CPML_regions
  logical, dimension(:), allocatable :: is_CPML

  ! moho surface (optional)
  integer :: nspec2D_moho
  integer, dimension(:), allocatable :: ibelm_moho
  integer, dimension(:,:), allocatable :: nodes_ibelm_moho

  character(len=256)  :: prname

  logical, dimension(:), allocatable :: mask_nodes_elmnts
  integer, dimension(:), allocatable :: used_nodes_elmnts

#ifdef USE_METIS_INSTEAD_OF_SCOTCH
  integer :: edgecut,wgtflag,numflag
  integer, dimension(0:4) :: dummy_options
  integer, dimension(1) :: dummy_array
#else
  double precision, dimension(SCOTCH_GRAPHDIM)  :: scotchgraph
  double precision, dimension(SCOTCH_STRATDIM)  :: scotchstrat
!!!!!! character(len=256), parameter :: scotch_strategy='b{job=t,map=t,poli=S,sep=h{pass=30}}'
#endif
  integer  :: ier,idummy

  !pll
  double precision , dimension(:,:), allocatable :: mat_prop
  integer :: count_def_mat,count_undef_mat,imat
  character (len=30), dimension(:,:), allocatable :: undef_mat_prop

! default mesh file directory
  character(len=256) :: localpath_name
  character(len=256) :: outputpath_name

  integer :: aniso_flag,idomain_id
  double precision :: vp,vs,rho,qkappa,qmu
! poroelastic parameters read in a new file
  double precision :: rhos,rhof,phi,tort,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,eta,mufr

  integer, parameter :: IIN_database = 15

  double precision :: DT
  double precision :: HDUR_MOVIE,OLSEN_ATTENUATION_RATIO,f0_FOR_PML
  integer :: NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP, &
            UTM_PROJECTION_ZONE,SIMULATION_TYPE,NGNOD,NGNOD2D
  integer :: NSOURCES,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY
  integer :: NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO,MOVIE_TYPE
  logical :: MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
            USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION
  logical :: ATTENUATION,USE_OLSEN_ATTENUATION,PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE, &
            APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,USE_FORCE_POINT_SOURCE,FULL_ATTENUATION_SOLID
  logical :: STACEY_ABSORBING_CONDITIONS,SAVE_FORWARD,STACEY_INSTEAD_OF_FREE_SURFACE
  logical :: ANISOTROPY,SAVE_MESH_FILES,USE_RICKER_TIME_FUNCTION,PRINT_SOURCE_TIME_FUNCTION
  character(len=256) LOCAL_PATH,TOMOGRAPHY_PATH,TRAC_PATH
  integer :: IMODEL

  contains

  !----------------------------------------------------------------------------------------------
  ! reads in mesh files
  !----------------------------------------------------------------------------------------------
  subroutine read_mesh_files

    implicit none

    character(len=256)  :: line
    logical :: use_poroelastic_file
    integer(long) :: nspec_long
    integer :: inode

  ! reads node coordinates
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/nodes_coords_file',&
          status='old', form='formatted', iostat = ier)
    if( ier /= 0 ) then
      print*,'could not open file:',localpath_name(1:len_trim(localpath_name))//'/nodes_coords_file'
      stop 'error file open'
    endif
    read(98,*) nnodes
    if( nnodes < 1 ) stop 'error: nnodes < 1'
    allocate(nodes_coords(3,nnodes),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_coords'
    do inode = 1, nnodes
    ! format: #id_node #x_coordinate #y_coordinate #z_coordinate
      read(98,*) num_node, nodes_coords(1,num_node), nodes_coords(2,num_node), nodes_coords(3,num_node)
    enddo
    close(98)
    print*, 'total number of nodes: '
    print*, '  nnodes = ', nnodes

    ! reads mesh elements indexing
    !(CUBIT calls this the connectivity, guess in the sense that it connects with the points index in
    ! the global coordinate file "nodes_coords_file"; it doesn't tell you which point is connected with others)
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/mesh_file', &
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening mesh_file'
    read(98,*) nspec_long

    ! debug check size limit
    if( nspec_long > 2147483646 ) then
      print *,'size exceeds integer 4-byte limit: ',nspec_long
      print*,'bit size fortran: ',bit_size(nspec)
      stop 'error number of elements too large'
    endif

    ! sets number of elements (integer 4-byte)
    nspec = nspec_long

    if( nspec < 1 ) stop 'error: nspec < 1'
    allocate(elmnts(NGNOD,nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array elmnts'
    do ispec = 1, nspec
      ! format: # element_id  #id_node1 ... #id_node8
      !      or # element_id  #id_node1 ... #id_node27

      ! note: be aware that here we can have different node ordering for a cube element;
      !          the ordering from Cubit files might not be consistent for multiple volumes, or uneven, unstructured grids
      !
      !          here our code assumes that element ordering is:
      !          at the bottom of the element, anticlock-wise, i.e.
      !             point 1 = (0,0,0), point 2 = (0,1,0), point 3 = (1,1,0), point 4 = (1,0,0)
      !          then top (positive z-direction) of element
      !             point 5 = (0,0,1), point 6 = (0,1,1), point 7 = (1,1,1), point 8 = (1,0,1)

      read(98,*,iostat=ier) num_elmnt,(elmnts(inode,num_elmnt), inode=1,NGNOD)

      if( ier /= 0 ) then
        print *,'error while attempting to read ',NGNOD,'element data values from the mesh file'
        if(NGNOD == 8) print *,'check if your mesh file is indeed composed of HEX8 elements'
        if(NGNOD == 27) print *,'check if your mesh file is indeed composed of HEX27 elements'
        stop 'error reading element data from the mesh file'
      endif

      if((num_elmnt > nspec) .or. (num_elmnt < 1) )  stop "ERROR : Invalid mesh file."

    enddo
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
      ! format: #id_element #flag
      ! note: be aware that elements may not be sorted in materials_file
      read(98,*) num_mat,mat(1,num_mat)
      if((num_mat > nspec) .or. (num_mat < 1) ) stop "ERROR : Invalid mat file."
    enddo
    close(98)

  ! reads material definitions
  !
  ! note: format of nummaterial_velocity_file must be
  !
  ! #(1)material_domain_id #(2)material_id  #(3)rho  #(4)vp   #(5)vs   #(6)Q_kappa   #(7)Q_mu  #(8)anisotropy_flag
  !
  ! where
  !     material_domain_id : 1=acoustic / 2=elastic / 3=poroelastic
  !     material_id        : number of material/volume
  !     rho                : density
  !     vp                 : P-velocity
  !     vs                 : S-velocity
  !     Q_kappa            : 9999 = no Q_kappa attenuation
  !     Q_mu               : 9999 = no Q_mu attenuation
  !     anisotropy_flag    : 0=no anisotropy/ 1,2,.. check with implementation in aniso_model.f90
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
       endif
       read(98,*,iostat=ier) idummy,num_mat
    enddo
    close(98)
    print*, '  defined = ',count_def_mat, 'undefined = ',count_undef_mat
    ! check with material flags
    if( count_def_mat > 0 .and. maxval(mat(1,:)) > count_def_mat ) then
      print*,'error material definitions:'
      print*,'  materials associated in materials_file:',maxval(mat(1,:))
      print*,'  larger than defined materials in nummaterial_velocity_file:',count_def_mat
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
       !
       ! reads lines until it reaches a defined material
       num_mat = -1
       do while( num_mat < 0 .and. ier == 0)
         read(98,'(A256)',iostat=ier) line
         read(line,*) idomain_id,num_mat
       enddo
       if( ier /= 0 ) stop 'error reading in defined materials in nummaterial_velocity_file'

       ! reads in defined material properties
       read(line,*) idomain_id,num_mat,rho,vp,vs,qkappa,qmu,aniso_flag

       ! sanity check: Q factor cannot be equal to zero, thus convert to 9999 to indicate no attenuation
       ! if users have used 0 to indicate that instead
       if(qkappa <= 0.000001) qkappa = 9999.
       if(qmu <= 0.000001) qmu = 9999.

       ! checks material_id bounds
       if(num_mat < 1 .or. num_mat > count_def_mat)  stop "ERROR : Invalid nummaterial_velocity_file file."

       if(idomain_id == 1 .or. idomain_id == 2) then ! material is elastic or acoustic

         ! check that the S-wave velocity is zero if the material is acoustic
         if(idomain_id == 1 .and. vs >= 0.0001) &
                stop 'acoustic material defined with a non-zero shear-wave velocity Vs, exiting...'

         ! check that the S-wave velocity is not zero if the material is elastic
         if(idomain_id == 2 .and. vs <= 0.0001) &
                stop '(visco)elastic material defined with a zero shear-wave velocity Vs, exiting...'

         mat_prop(1,num_mat) = rho
         mat_prop(2,num_mat) = vp
         mat_prop(3,num_mat) = vs
         mat_prop(4,num_mat) = qmu
         mat_prop(5,num_mat) = aniso_flag
         mat_prop(6,num_mat) = idomain_id
         mat_prop(7,num_mat) = qkappa  ! this one is not stored next to qmu for historical reasons, because it was added later

       else if(idomain_id == 3) then ! material is poroelastic

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

       else
         stop 'idomain_id must be 1, 2 or 3 for acoustic, elastic or poroelastic in nummaterial_velocity_file'

       endif ! of if(idomain_id == ...)

    enddo

    ! reads in undefined material properties
    rewind(98,iostat=ier) ! back to the beginning of the file
    do imat=1,count_undef_mat
       !  undefined materials: have to be listed in decreasing order of material_id (start with -1, -2, etc...)
       !  format:
       !   - for interfaces
       !    #(6) material_domain_id #(1) material_id(<0) #(2) type_name (="interface")
       !     #(3) material_id_for_material_below #(4) material_id_for_material_above
       !        example:     2 -1 interface 1 2
       !   - for tomography models
       !    #(6) material_domain_id #(1) material_id(<0) #(2) type_name (="tomography")
       !     #(3) block_name (="elastic") #(4) file_name
       !        example:     2 -1 tomography elastic tomography_model.xyz
       ! reads lines until it reaches a defined material
       num_mat = 1
       do while( num_mat >= 0 .and. ier == 0 )
         read(98,'(A256)',iostat=ier) line
         read(line,*) idomain_id,num_mat
       enddo
       if( ier /= 0 ) stop 'error reading in undefined materials in nummaterial_velocity_file'

       ! checks if interface or tomography definition
       read(line,*) undef_mat_prop(6,imat),undef_mat_prop(1,imat),undef_mat_prop(2,imat)
       read(undef_mat_prop(1,imat),*) num_mat
       if( trim(undef_mat_prop(2,imat)) == 'interface' ) then
         ! line will have 5 arguments, e.g.: 2 -1 interface 1 2
         read(line,*) undef_mat_prop(6,imat),undef_mat_prop(1,imat),undef_mat_prop(2,imat),&
                     undef_mat_prop(3,imat),undef_mat_prop(4,imat)
         undef_mat_prop(5,imat) = "0" ! dummy value
       else if( trim(undef_mat_prop(2,imat)) == 'tomography' ) then
         ! line will have 6 arguments, e.g.: 2 -1 tomography elastic tomography_model.xyz 1
         read(line,*) undef_mat_prop(6,imat),undef_mat_prop(1,imat),undef_mat_prop(2,imat),&
                        undef_mat_prop(3,imat),undef_mat_prop(4,imat)
         undef_mat_prop(5,imat) = "0" ! dummy value
       else
         stop "ERROR: invalid line in nummaterial_velocity_file for undefined material"
       endif

       ! checks material_id
       if( trim(undef_mat_prop(2,imat)) == 'interface' .or. trim(undef_mat_prop(2,imat)) == 'tomography' ) then
          if(num_mat > 0 .or. -num_mat > count_undef_mat)  &
               stop "ERROR : Invalid nummaterial_velocity_file for undefined materials."
          if(num_mat /= -imat)  &
               stop "ERROR : Invalid material_id in nummaterial_velocity_file for undefined materials."
       endif

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
    enddo
    if( use_poroelastic_file ) close(97)
    close(98)

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
            else if( trim(undef_mat_prop(2,imat)) == 'tomography' ) then
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
! if the file does not exist then define the number of Stacey elements as zero for this face;
! beware that these files can also be used to set Dirichlet boundary conditions on the outer edges of CPML
! absorbing layers for elastic elements, not only for Stacey; thus these files may exist and be non-empty
! even when STACEY_ABSORBING_CONDITIONS is false
    if( ier /= 0 ) then
      nspec2D_xmin = 0
    else
      read(98,*) nspec2D_xmin
    endif
! 33333333333333333333333333333333333333333
! an array of size 0 is a valid object in Fortran 90, i.e. the array is then considered as allocated
! and can thus for instance be used as an argument in a call to a subroutine without giving any error
! even when full range and pointer checking is used in the compiler options;
! thus here the idea is that if some of the absorbing files do not exist because there are no absorbing
! conditions for this mesh then the array is created nonetheless, but with a dummy size of 0
    allocate(ibelm_xmin(nspec2D_xmin),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_xmin'
    allocate(nodes_ibelm_xmin(NGNOD2D,nspec2D_xmin),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_xmin'
    do ispec2D = 1,nspec2D_xmin
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      ! note: ordering for CUBIT seems such that the normal of the face points outward of the element the face belongs to;
      !         in other words, nodes are in increasing order such that when looking from within the element outwards,
      !         they are ordered clockwise
      !
      !          doesn't necessarily have to start on top-rear, then bottom-rear, bottom-front, and finally top-front i.e.:
      !          point 1 = (0,1,1), point 2 = (0,1,0), point 3 = (0,0,0), point 4 = (0,0,1)
      read(98,*) ibelm_xmin(ispec2D), (nodes_ibelm_xmin(inode,ispec2D), inode=1,NGNOD2D)
    enddo
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
    allocate(nodes_ibelm_xmax(NGNOD2D,nspec2D_xmax),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_xmax'
    do ispec2D = 1,nspec2D_xmax
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_xmax(ispec2D), (nodes_ibelm_xmax(inode,ispec2D), inode=1,NGNOD2D)
    enddo
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
    allocate(nodes_ibelm_ymin(NGNOD2D,nspec2D_ymin),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_ymin'
    do ispec2D = 1,nspec2D_ymin
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_ymin(ispec2D), (nodes_ibelm_ymin(inode,ispec2D), inode=1,NGNOD2D)
    enddo
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
    allocate(nodes_ibelm_ymax(NGNOD2D,nspec2D_ymax),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_ymax'
    do ispec2D = 1,nspec2D_ymax
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_ymax(ispec2D), (nodes_ibelm_ymax(inode,ispec2D), inode=1,NGNOD2D)
    enddo
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
    allocate(nodes_ibelm_bottom(NGNOD2D,nspec2D_bottom),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_bottom'
    do ispec2D = 1,nspec2D_bottom
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_bottom(ispec2D), (nodes_ibelm_bottom(inode,ispec2D), inode=1,NGNOD2D)
    enddo
    close(98)
    print*, '  nspec2D_bottom = ', nspec2D_bottom

  ! reads in free_surface boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/free_or_absorbing_surface_file_zmax', &
          status='old', form='formatted',iostat=ier)
    if( ier /= 0 ) then
      nspec2D_top = 0
    else
      read(98,*) nspec2D_top
    endif
    allocate(ibelm_top(nspec2D_top),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_top'
    allocate(nodes_ibelm_top(NGNOD2D,nspec2D_top),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_top'
    do ispec2D = 1,nspec2D_top
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_top(ispec2D), (nodes_ibelm_top(inode,ispec2D), inode=1,NGNOD2D)
    enddo
    close(98)
    print*, '  nspec2D_top = ', nspec2D_top

! 33333333333333333333333333333333333333333
! an array of size 0 is a valid object in Fortran 90, i.e. the array is then considered as allocated
! and can thus for instance be used as an argument in a call to a subroutine without giving any error
! even when full range and pointer checking is used in the compiler options;
! thus here the idea is that if some of the absorbing files do not exist because there are no absorbing
! conditions for this mesh then the array is created nonetheless, but with a dummy size of 0

  ! reads in absorbing_cpml boundary file
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_cpml_file', &
         status='old', form='formatted',iostat=ier)
! if the file does not exist but if there are PML_CONDITIONS then stop
    if( ier /= 0 .and. PML_CONDITIONS) &
        stop 'error: PML_CONDITIONS is set to true but file absorbing_cpml_file does not exist'
! if the file does not exist or if there are no PML_CONDITIONS then define the number of CPML elements as zero
    if( ier /= 0 .or. .not. PML_CONDITIONS) then
       nspec_cpml = 0
    else
       read(98,*) nspec_cpml
    endif

! sanity check
    if( PML_CONDITIONS .and. nspec_cpml <= 0) &
        stop 'error: PML_CONDITIONS is set to true but nspec_cpml <= 0 in file absorbing_cpml_file'

    ! C-PML spectral elements global indexing
    allocate(CPML_to_spec(nspec_cpml),stat=ier)
    if(ier /= 0) stop 'error allocating array CPML_to_spec'
    ! C-PML regions (see below)
    allocate(CPML_regions(nspec_cpml),stat=ier)
    if(ier /= 0) stop 'error allocating array CPML_regions'
    do ispec_CPML=1,nspec_cpml
       ! elements are stored with #id_cpml_regions increasing order:
       !
       ! #id_cpml_regions = 1 : X_surface C-PML
       ! #id_cpml_regions = 2 : Y_surface C-PML
       ! #id_cpml_regions = 3 : Z_surface C-PML
       ! #id_cpml_regions = 4 : XY_edge C-PML
       ! #id_cpml_regions = 5 : XZ_edge C-PML
       ! #id_cpml_regions = 6 : YZ_edge C-PML
       ! #id_cpml_regions = 7 : XYZ_corner C-PML
       !
       ! format: #id_cpml_element #id_cpml_regions
       read(98,*) CPML_to_spec(ispec_CPML), CPML_regions(ispec_CPML)
    enddo
    close(98)
    if( nspec_cpml > 0 ) print*, '  nspec_cpml = ', nspec_cpml

    ! sets mask of C-PML elements for all elements in this partition
    allocate(is_CPML(nspec),stat=ier)
    if(ier /= 0) stop 'error allocating array is_CPML'
    is_CPML(:) = .false.
    do ispec_CPML=1,nspec_cpml
       if( (CPML_regions(ispec_CPML)>=1) .and. (CPML_regions(ispec_CPML)<=7) ) then
          is_CPML(CPML_to_spec(ispec_CPML)) = .true.
       endif
    enddo

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
    allocate(nodes_ibelm_moho(NGNOD2D,nspec2D_moho),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_ibelm_moho'
    do ispec2D = 1,nspec2D_moho
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_moho(ispec2D), (nodes_ibelm_moho(inode,ispec2D), inode=1,NGNOD2D)
    enddo
    close(98)
    if( nspec2D_moho > 0 ) print*, '  nspec2D_moho = ', nspec2D_moho

    call read_fault_files(localpath_name)
    if (ANY_FAULT) then
      call save_nodes_coords(nodes_coords,nnodes)
      call close_faults(nodes_coords,nnodes)
    endif

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
      do inode = 1, NGNOD
        mask_nodes_elmnts(elmnts(inode,ispec)) = .true.
        used_nodes_elmnts(elmnts(inode,ispec)) = used_nodes_elmnts(elmnts(inode,ispec)) + 1
      enddo
    enddo
    print *, 'node valence:'
    print *, '  min = ',minval(used_nodes_elmnts(:)),' max = ', maxval(used_nodes_elmnts(:))
    do inode = 1, nnodes
      if (.not. mask_nodes_elmnts(inode)) then
        stop 'ERROR: found some unused nodes (weird, but not necessarily fatal; your mesher may have created extra nodes).'
      endif
    enddo

! max number of elements that contain the same node
    nsize = maxval(used_nodes_elmnts(:))

! majoration (overestimate) of the maximum number of neighbours per element
!! DK DK nfaces is a constant equal to 6 (number of faces of a cube).
!! DK DK I have no idea how this formula works; it was designed by Nicolas Le Goff
    sup_neighbour = NGNOD_EIGHT_CORNERS * nsize - (NGNOD_EIGHT_CORNERS + (NGNOD_EIGHT_CORNERS/2 - 1)*nfaces)

    ! debug check size limit
!! DK DK this check will likely fail because sup_neighbour itself may become negative if going over the 4-byte integer limit;
!! DK DK but this should never happen in practice (by far)...
    if( sup_neighbour > 2147483646 ) then
      print *,'size exceeds integer 4-byte limit: ',sup_neighbour,nsize
      print *,'bit size fortran: ',bit_size(sup_neighbour)
      stop 'ERROR: sup_neighbour is too large'
    endif

    print *, '  nsize = ',nsize, 'sup_neighbour = ', sup_neighbour

  end subroutine check_valence

  !----------------------------------------------------------------------------------------------
  ! divides model into partitions using scotch library functions
  !----------------------------------------------------------------------------------------------

  subroutine scotch_partitioning

    implicit none
    ! local parameters
    integer, dimension(:), allocatable  :: num_material
    integer :: ier

    ! starts from 0
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DK DK added this in Oct 2012 to see if we first do 4 and then 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_TWO_CALLS_TO_mesh2dual
    ncommonnodes = NGNOD2D_FOUR_CORNERS
#else
    ncommonnodes = 1
#endif
    call mesh2dual_ncommonnodes(nspec, nnodes, nsize, sup_neighbour, elmnts, xadj, adjncy, nnodes_elmnts, &
         nodes_elmnts, max_neighbour, ncommonnodes, NGNOD)

#ifdef USE_TWO_CALLS_TO_mesh2dual
    print*, 'mesh2dual first call:'
#else
    print*, 'mesh2dual:'
#endif
    print*, '  max_neighbour = ',max_neighbour

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DK DK added this in Oct 2012 to see if we first do 4 and then 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! DK DK Oct 2012: added this safety test
    if(max_neighbour > sup_neighbour) stop 'found max_neighbour > sup_neighbour in domain decomposition'

    nb_edges = xadj(nspec+1)

    ! allocates & initializes partioning of elements
    allocate(part(1:nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array part'
    part(:) = -1

    ! initializes
    ! elements load array
    allocate(elmnts_load(1:nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array elmnts_load'

    ! gets materials id associations
    allocate(num_material(1:nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array num_material'
    ! note: num_material can be negative for tomographic material elements
    !       (which are counted then as elastic elements)
    num_material(:) = mat(1,:)

!! DK DK Oct 2012: this should include CPML weights as well in the future
    ! uniform load by default
    elmnts_load(:) = ACOUSTIC_LOAD
    ! then in case of acoustic/elastic/poro simulation, assign different weights to elements accordingly
    call acoustic_elastic_poro_load(elmnts_load,nspec,count_def_mat,count_undef_mat, &
                                  num_material,mat_prop,undef_mat_prop,ATTENUATION)

#ifdef USE_METIS_INSTEAD_OF_SCOTCH

    ! METIS partitioning

    edgecut = 0
    dummy_options(0:4) = 0
    wgtflag = 2
    numflag = 0
!! DK DK Metis4.0 syntax
!! DK DK we should find a way of testing the return code of this function in Fortran
    call METIS_PartGraphKway(nspec, xadj, adjncy, elmnts_load, dummy_array, wgtflag, numflag, nparts, &
                                  dummy_options, edgecut, part)
    print *
    print *,'total edgecut created by METIS = ',edgecut
    print *

#else

    ! SCOTCH partitioning

    ! we use default strategy for partitioning, thus omit specifing explicit strategy.

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
    !!This will make you independent from further changes
    !!(improvements) in the strategy syntax.
    !!And you should see an improvement in performance, too,
    !!as your hand-made strategy did not make use of the
    !!multi-level framework.

    call scotchfstratinit (scotchstrat(1), ier)
     if (ier /= 0) then
       stop 'ERROR : MAIN : Cannot initialize strategy'
    endif

    !call scotchfstratgraphmap (scotchstrat(1), trim(scotch_strategy), ier)
    ! if (ier /= 0) then
    !   stop 'ERROR : MAIN : Cannot build strategy'
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
       stop 'ERROR : MAIN : Cannot destroy strategy'
    endif

#endif

    ! re-partitioning puts poroelastic-elastic coupled elements into same partition
    !  integer  :: nfaces_coupled
    !  integer, dimension(:,:), pointer  :: faces_coupled
    ! TODO: supposed to rebalance, but currently broken
!! DK DK added this because poroelastic repartitioning routine of Christina Morency is currently broken
! implement mesh repartitioning of poroelastic-elastic interface
! (the risk being to break the nice load balancing created by the domain decomposer for high-performance computing)
    if(PORO_INTERFACE_REPARTITIONING) &
      call poro_elastic_repartitioning (nspec, nnodes, elmnts, &
                       count_def_mat, num_material , mat_prop, &
                       sup_neighbour, nsize, &
                       nparts, part, NGNOD)

    deallocate(num_material)

  ! re-partitioning transfers two coupled elements on fault side 1 and side 2 to the same partition
    if (ANY_FAULT) call fault_repartition (nspec, nnodes, elmnts, nsize, nparts, part, NGNOD, nodes_coords)


    ! re-partitioning puts moho-surface coupled elements into same partition
! (the risk being to break the nice load balancing created by the domain decomposer for high-performance computing)
    if(SAVE_MOHO_MESH) call moho_surface_repartitioning (nspec, nnodes, elmnts, &
                                        sup_neighbour, nsize, nparts, part, &
                                        nspec2D_moho,ibelm_moho,nodes_ibelm_moho, NGNOD, NGNOD2D)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DK DK added this in Oct 2012 to see if we first do 4 and then 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_TWO_CALLS_TO_mesh2dual

    ncommonnodes = 1
    call mesh2dual_ncommonnodes(nspec, nnodes, nsize, sup_neighbour, elmnts, xadj, adjncy, nnodes_elmnts, &
         nodes_elmnts, max_neighbour, ncommonnodes, NGNOD)

    print*, 'mesh2dual second call:'
    print*, '  max_neighbour = ',max_neighbour

!! DK DK Oct 2012: added this safety test
    if(max_neighbour > sup_neighbour) stop 'found max_neighbour > sup_neighbour in domain decomposition'

    nb_edges = xadj(nspec+1)

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DK DK added this in Oct 2012 to see if we first do 4 and then 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
                             nparts, NGNOD)

    ! obsolete: from when we wanted acoustic/elastic boundaries NOT to be separated into different MPI partitions
    ! call build_interfaces_no_ac_el_sep(nspec, sup_neighbour, part, elmnts, &
    !                           xadj, adjncy, tab_interfaces, &
    !                           tab_size_interfaces, ninterfaces, &
    !                           count_def_mat, mat_prop(3,:), mat(1,:), nparts, NGNOD)

  end subroutine scotch_partitioning

  !----------------------------------------------------------------------------------------------
  ! writes out new Databases files for each partition
  !----------------------------------------------------------------------------------------------

  subroutine write_mesh_databases

    implicit none

    integer :: ier

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

       call write_partition_database(IIN_database, ipart, nspec_local, nspec, elmnts, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, mat, NGNOD, 1)

       !debug
       !print*, ipart,": nspec_local=",nspec_local, " nnodes_local=", nnodes_loc

       ! writes out node coordinate locations
       write(IIN_database) nnodes_loc

       call write_glob2loc_nodes_database(IIN_database, ipart, nnodes_loc, nodes_coords,&
                                  glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, nnodes, 2)

       call write_material_props_database(IIN_database,count_def_mat,count_undef_mat, &
                                  mat_prop, undef_mat_prop)

       ! writes out spectral element indices
       write(IIN_database) nspec_local
       call write_partition_database(IIN_database, ipart, nspec_local, nspec, elmnts, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, mat, NGNOD, 2)

       ! writes out absorbing/free-surface boundaries
       call write_boundaries_database(IIN_database, ipart, nspec, nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, &
                                  nspec2D_ymax, nspec2D_bottom, nspec2D_top, &
                                  ibelm_xmin, ibelm_xmax, ibelm_ymin, &
                                  ibelm_ymax, ibelm_bottom, ibelm_top, &
                                  nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
                                  nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, NGNOD2D)

       ! writes out C-PML elements indices, CPML-regions and thickness of C-PML layer
       call write_cpml_database(IIN_database, ipart, nspec, nspec_cpml, CPML_to_spec, &
            CPML_regions, is_CPML, glob2loc_elmnts, part)

       ! gets number of MPI interfaces
       call Write_interfaces_database(IIN_database, tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
                                  my_ninterface, my_interfaces, my_nb_interfaces, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, 1, nparts)

       ! writes out MPI interfaces elements
       !print*,' my interfaces:',my_ninterface,maxval(my_nb_interfaces)
       if( my_ninterface == 0 ) then
        write(IIN_database) my_ninterface, 0       ! avoids problem with maxval for empty array my_nb_interfaces
       else
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
                                  nspec2D_moho, ibelm_moho, nodes_ibelm_moho, NGNOD2D)

       close(IIN_database)


       ! write fault database
       if (ANY_FAULT) then
          write(prname, "(i6.6,'_Database_fault')") ipart
          open(unit=16,file=outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname,&
               status='replace', action='write', form='unformatted', iostat = ier)
          if( ier /= 0 ) then
            print*,'error file open:',outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname
            print*
            print*,'check if path exists:',outputpath_name(1:len_trim(outputpath_name))
            stop
          endif
          call write_fault_database(16, ipart, nspec, &
                                    glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                    glob2loc_nodes, part)
          !write(16,*) nnodes_loc
          write(16) nnodes_loc
          call write_glob2loc_nodes_database(16, ipart, nnodes_loc, nodes_coords_open,&
                                  glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, nnodes, 2)
          close(16)
       endif


    enddo

    ! cleanup
    deallocate(CPML_to_spec,stat=ier); if( ier /= 0 ) stop 'error deallocating array CPML_to_spec'
    deallocate(CPML_regions,stat=ier); if( ier /= 0 ) stop 'error deallocating array CPML_regions'
    deallocate(is_CPML,stat=ier); if( ier /= 0 ) stop 'error deallocating array is_CPML'

    print*, 'partitions: '
    print*, '  num = ',nparts
    print*
    print*, 'Databases files in directory: ',outputpath_name(1:len_trim(outputpath_name))
    print*, 'finished successfully'
    print*

  end subroutine write_mesh_databases

end module decompose_mesh

