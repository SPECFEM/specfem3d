!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

module module_mesh

  use constants, only: NDIM,MAX_STRING_LEN,IIN_DB,IIN_DB2
  use shared_parameters, only: NGNOD,NGNOD2D,ATTENUATION,PML_CONDITIONS
  use fault_scotch

  ! elements
  integer                                                    :: nspec, nspec_glob
  integer,               dimension(:,:),        allocatable  :: elmnts, elmnts_glob
  integer,               dimension(:,:),        allocatable  :: mat

  integer                                                    :: nspec_part_boundaries
  integer,               dimension(:,:),        allocatable  :: elmnts_part_boundaries
  integer,               dimension(:),          allocatable  :: iboundary

  ! vertices
  integer                                                    :: nnodes, nnodes_glob
  double precision,      dimension(:,:),         allocatable :: nodes_coords, nodes_coords_glob
  double precision,      dimension(:,:),         allocatable :: nodes_coords_open_loc

  ! boundaries
  integer                                                    :: ispec2D
  integer                                                    :: nspec2D_xmin, nspec2D_xmax
  integer                                                    :: nspec2D_ymin, nspec2D_ymax
  integer                                                    :: nspec2D_bottom, nspec2D_top
  integer,               dimension(:),           allocatable :: ibelm_xmin, ibelm_xmax, ibelm_ymin
  integer,               dimension(:),           allocatable :: ibelm_ymax, ibelm_bottom, ibelm_top
  integer,               dimension(:,:),         allocatable :: nodes_ibelm_xmin, nodes_ibelm_xmax
  integer,               dimension(:,:),         allocatable :: nodes_ibelm_ymin, nodes_ibelm_ymax
  integer,               dimension(:,:),         allocatable :: nodes_ibelm_bottom, nodes_ibelm_top

  ! C-PML absorbing boundary conditions
  integer                                                    :: ispec_cpml
  integer                                                    :: nspec_cpml
  integer,               dimension(:),           allocatable :: cpml_to_spec, cpml_regions
  logical,               dimension(:),           allocatable :: is_CPML

  ! moho surface (optional)
  integer                                                    :: nspec2D_moho
  integer,              dimension(:),            allocatable :: ibelm_moho
  integer,              dimension(:,:),          allocatable :: nodes_ibelm_moho


  ! material properties
  integer                                                    :: count_def_mat,count_undef_mat,imat
  integer,               dimension(:),           allocatable :: num_material
  double precision ,     dimension(:,:),         allocatable :: mat_prop
  character(len=MAX_STRING_LEN), dimension(:,:), allocatable :: undef_mat_prop

  ! elements load
  double precision,      dimension(:),           allocatable :: load_elmnts

  !! acoustic-elastic-poroelastic as well as CPML load balancing:
  !! we define here the relative cost of all types of spectral elements used in the code.
  double precision, parameter :: ACOUSTIC_LOAD     = 1.d0
  double precision, parameter :: ELASTIC_LOAD      = 4.1d0
  double precision, parameter :: VISCOELASTIC_LOAD = 5.9d0
  double precision, parameter :: POROELASTIC_LOAD  = 8.1d0


  ! default mesh file directory
  character(len=MAX_STRING_LEN) :: localpath_name

  ! useful kind types for short integer (4 bytes) and long integers (8 bytes)
  integer, parameter :: short = 4, long = 8

contains

  !----------------------------------------------------------------------------------------------
  ! reads in mesh files
  !----------------------------------------------------------------------------------------------
  subroutine read_mesh_files()

  ! note: this routine will be called by the parallel decompose mesher.
  !       most of the routine is a duplicate from the routine read_mesh_files() in src/decompose_mesh/decompose_mesh.F90

    implicit none

    character(len=MAX_STRING_LEN)                            :: line
    logical                                                  :: use_poroelastic_file
    integer(long)                                            :: nspec_long
    integer                                                  :: inode, num_node, ispec, num_elmnt
    integer                                                  :: num_mat
    integer                                                  :: ier, idummy
    integer                                                  :: aniso_flag,idomain_id
    double precision                                         :: vp,vs,rho,qkappa,qmu
    ! poroelastic parameters read in a new file
    double precision :: rhos,rhof,phi,tort,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,eta,mufr
    ! position
    double precision :: x_coord,y_coord,z_coord
    integer, dimension(NGNOD) :: elmnts_ids
    integer :: mat_flag

    localpath_name='./MESH'

    ! reads node coordinates
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/nodes_coords_file', &
         status='old', form='formatted', iostat = ier)
    if (ier /= 0) then
       print *,'could not open file:',localpath_name(1:len_trim(localpath_name))//'/nodes_coords_file'
       stop 'Error opening file nodes_coords_file'
    endif
    read(IIN_DB,*) nnodes_glob

    if (nnodes_glob < 1) stop 'Error: nnodes_glob < 1'

    allocate(nodes_coords_glob(NDIM,nnodes_glob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 117')
    if (ier /= 0) stop 'Error allocating array nodes_coords'
    nodes_coords_glob(:,:) = 0.d0

    do inode = 1, nnodes_glob
      ! format: #id_node #x_coordinate #y_coordinate #z_coordinate
      read(IIN_DB,*) num_node,x_coord,y_coord,z_coord

      ! stores coordinates
      nodes_coords_glob(1,num_node) = x_coord
      nodes_coords_glob(2,num_node) = y_coord
      nodes_coords_glob(3,num_node) = z_coord

      ! for parallel meshing
      if (mod(inode,100000) == 0) then
        write(27,'(2i10)') num_node/100000, nnodes_glob/100000
      endif
    enddo
    close(IIN_DB)
    write(27,*) 'total number of nodes: '
    write(27,*) '  nnodes = ', nnodes_glob

    ! reads mesh elements indexing
    !(CUBIT calls this the connectivity, guess in the sense that it connects with the points index in
    ! the global coordinate file "nodes_coords_file"; it doesn't tell you which point is connected with others)
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/mesh_file', &
         status='old', form='formatted',iostat=ier)
    if (ier /= 0) stop 'Error opening mesh_file'
    read(IIN_DB,*) nspec_long

    ! debug check size limit
    if (nspec_long > 2147483646) then
       print *,'size exceeds integer 4-byte limit: ',nspec_long
       print *,'bit size Fortran: ',bit_size(nspec)
       stop 'Error number of elements too large'
    endif

    ! sets number of elements (integer 4-byte)
    nspec_glob = nspec_long

    if (nspec_glob < 1) stop 'Error: nspec_glob < 1'

    allocate(elmnts_glob(NGNOD,nspec_glob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 118')
    if (ier /= 0) stop 'Error allocating array elmnts'
    elmnts_glob(:,:) = 0

    do ispec = 1, nspec_glob
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

      read(IIN_DB,*,iostat=ier) num_elmnt,(elmnts_ids(inode), inode=1,NGNOD)

      ! stores node ids
      elmnts_glob(:,num_elmnt) = elmnts_ids(:)

       if (ier /= 0) then
          print *,'Error while attempting to read ',NGNOD,'element data values from the mesh file'
          if (NGNOD == 8) print *,'check if your mesh file is indeed composed of HEX8 elements'
          if (NGNOD == 27) print *,'check if your mesh file is indeed composed of HEX27 elements'
          stop 'Error reading element data from the mesh file'
       endif

       if ((num_elmnt > nspec_glob) .or. (num_elmnt < 1))  stop "Error : Invalid mesh_file"

       ! for parallel meshing
       if (mod(ispec,100000) == 0) then
          write(27,'(2i10)') ispec/100000, nspec_glob/100000
       endif
    enddo
    close(IIN_DB)
    write(27,*) 'total number of spectral elements:'
    write(27,*) '  nspec = ', nspec_glob

    ! reads material associations
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/materials_file', &
         status='old', form='formatted',iostat=ier)
    if (ier /= 0) stop 'Error opening materials_file'

    allocate(mat(2,nspec_glob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 119')
    if (ier /= 0) stop 'Error allocating array mat'
    mat(:,:) = 0

    do ispec = 1, nspec_glob
      ! format: #id_element #flag
      ! note: be aware that elements may not be sorted in materials_file
      read(IIN_DB,*) num_mat,mat_flag

      ! store material number
      mat(1,num_mat) = mat_flag

      if ((num_mat > nspec_glob) .or. (num_mat < 1)) stop "Error : Invalid materials_file"
    enddo
    close(IIN_DB)

    ! gets materials id associations
    allocate(num_material(1:nspec_glob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 120')
    if (ier /= 0) stop 'Error allocating array num_material'
    ! note: num_material can be negative for tomographic material elements
    !       (which are counted then as elastic elements)
    num_material(:) = mat(1,:)

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
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/nummaterial_velocity_file', &
         status='old', form='formatted',iostat=ier)
    if (ier /= 0) stop 'Error opening nummaterial_velocity_file'

    ! counts materials (defined/undefined)
    write(27,*) 'materials:'
    ier = 0
    do while (ier == 0)
       ! note: format #material_domain_id #material_id #...
       read(IIN_DB,'(A)',iostat=ier) line
       if (ier /= 0) exit

       ! skip empty/comment lines
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#' .or. line(1:1) == '!') cycle

       read(line,*,iostat=ier) idummy,num_mat
       if (ier /= 0) exit

       write(27,*) '  num_mat = ',num_mat
       ! checks non-zero material id
       if (num_mat == 0) stop 'Error in nummaterial_velocity_file: material id 0 found. Material ids must be non-zero.'
       if (num_mat > 0) then
          ! positive materials_id: velocity values will be defined
          count_def_mat = count_def_mat + 1
       else
          ! negative materials_id: undefined material properties yet
          count_undef_mat = count_undef_mat + 1
       endif
    enddo
    write(27,*) '  defined = ',count_def_mat, 'undefined = ',count_undef_mat

    ! check with material flags
    ! defined materials
    if (count_def_mat > 0 .and. maxval(mat(1,:)) > count_def_mat) then
       print *,'Error material definitions:'
       print *,'  materials associated in materials_file:',maxval(mat(1,:))
       print *,'  larger than defined materials in nummaterial_velocity_file:',count_def_mat
       stop 'Error positive material id exceeds bounds for defined materials'
    endif
    allocate(mat_prop(17,count_def_mat),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 121')
    if (ier /= 0) stop 'Error allocating array mat_prop'
    mat_prop(:,:) = 0.d0

    ! undefined materials
    if (count_undef_mat > 0 .and. minval(mat(1,:)) < -count_undef_mat) then
       print *,'Error material definitions:'
       print *,'  materials associated in materials_file:',minval(mat(1,:))
       print *,'  larger than undefined materials in nummaterial_velocity_file:',count_undef_mat
       stop 'Error negative material id exceeds bounds for undefined materials'
    endif
    allocate(undef_mat_prop(6,count_undef_mat),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 122')
    if (ier /= 0) stop 'Error allocating array undef_mat_prop'
    undef_mat_prop(:,:) = ''

    ! reads in defined material properties
    rewind(IIN_DB,iostat=ier)
    if (ier /= 0) stop 'Error rewinding nummaterial_velocity_file'

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
    open(unit=IIN_DB2, file=localpath_name(1:len_trim(localpath_name))//'/nummaterial_poroelastic_file', &
         status='old', form='formatted', iostat=ier)
    ! checks if we can use file
    if (ier /= 0) then
       use_poroelastic_file = .false.
       !stop 'Error opening nummaterial_poroelastic_file'
    else
       use_poroelastic_file = .true.
       write(27,*) '  poroelastic material file found'
    endif
    ier = 0

    ! note: entries in nummaterial_velocity_file can be an unsorted list of all
    !          defined materials (material_id > 0) and undefined materials (material_id < 0)
    do imat=1,count_def_mat
       ! material definitions
       !
       ! format: note that we save the arguments in a slightly different order in mat_prop(:,:)
       !
       ! reads lines until it reaches a defined material
       num_mat = -1
       do while (num_mat < 0 .and. ier == 0)
          read(IIN_DB,'(A)',iostat=ier) line
          if (ier /= 0) exit

          ! skip empty/comment lines
          if (len_trim(line) == 0) cycle
          if (line(1:1) == '#' .or. line(1:1) == '!') cycle

          read(line,*) idomain_id,num_mat
       enddo
       if (ier /= 0) stop 'Error reading in defined materials in nummaterial_velocity_file'

       ! reads in defined material properties
       read(line,*) idomain_id,num_mat,rho,vp,vs,qkappa,qmu,aniso_flag

       ! sanity check: Q factor cannot be equal to zero, thus convert to 9999 to indicate no attenuation
       ! if users have used 0 to indicate that instead
       if (qkappa <= 0.000001) qkappa = 9999.
       if (qmu <= 0.000001) qmu = 9999.

       ! checks material_id bounds
       if (num_mat < 1 .or. num_mat > count_def_mat) &
            stop "Error in nummaterial_velocity_file: material id invalid for defined materials."

       if (idomain_id == 1 .or. idomain_id == 2) then
          ! material is elastic or acoustic

          ! check that the S-wave velocity is zero if the material is acoustic
          if (idomain_id == 1 .and. vs >= 0.0001) &
               stop 'Error in nummaterial_velocity_file: acoustic material defined with a non-zero shear-wave velocity Vs.'

          ! check that the S-wave velocity is not zero if the material is elastic
          if (idomain_id == 2 .and. vs <= 0.0001) &
               stop 'Error in nummaterial_velocity_file: (visco)elastic material defined with a zero shear-wave velocity Vs.'

          mat_prop(1,num_mat) = rho
          mat_prop(2,num_mat) = vp
          mat_prop(3,num_mat) = vs
          mat_prop(4,num_mat) = qkappa
          mat_prop(5,num_mat) = qmu
          mat_prop(6,num_mat) = aniso_flag
          mat_prop(7,num_mat) = idomain_id

       else if (idomain_id == 3) then
          ! material is poroelastic

          if (use_poroelastic_file .eqv. .false.) &
               stop 'Error in nummaterial_velocity_file: poroelastic material requires nummaterial_poroelastic_file'

          ! reads poroelastic file line
          read(IIN_DB2,'(A)',iostat=ier) line
          if (ier /= 0) stop 'Error reading line in nummaterial_poroelastic_file'

          read(line,*) rhos,rhof,phi,tort,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,eta,mufr

          mat_prop(1,num_mat) = rhos
          mat_prop(2,num_mat) = rhof
          mat_prop(3,num_mat) = phi
          mat_prop(4,num_mat) = tort
          mat_prop(5,num_mat) = eta
          mat_prop(6,num_mat) = 0     ! aniso_flag
          mat_prop(7,num_mat) = idomain_id
          mat_prop(8,num_mat) = kxx
          mat_prop(9,num_mat) = kxy
          mat_prop(10,num_mat) = kxz
          mat_prop(11,num_mat) = kyy
          mat_prop(12,num_mat) = kyz
          mat_prop(13,num_mat) = kzz
          mat_prop(14,num_mat) = kappas
          mat_prop(15,num_mat) = kappaf
          mat_prop(16,num_mat) = kappafr
          mat_prop(17,num_mat) = mufr

       else
          stop 'Error in nummaterial_velocity_file: idomain_id must be 1, 2 or 3 for acoustic, elastic or poroelastic domains'

       endif ! of if (idomain_id == ...)

    enddo

    ! back to the beginning of the file
    rewind(IIN_DB,iostat=ier)
    if (ier /= 0) stop 'Error rewinding nummaterial_velocity_file'

    ! reads in undefined material properties
    do imat = 1,count_undef_mat
       !  undefined materials: have to be listed in decreasing order of material_id (start with -1, -2, etc...)
       !  format:
       !   - for interfaces
       !    #(6) material_domain_id #(1) material_id( < 0) #(2) type_name (="interface")
       !     #(3) material_id_for_material_below #(4) material_id_for_material_above
       !        example:     2 -1 interface 1 2
       !   - for tomography models
       !    #(6) material_domain_id #(1) material_id( < 0) #(2) type_name (="tomography")
       !     #(3) block_name (="elastic") #(4) file_name
       !        example:     2 -1 tomography elastic tomography_model.xyz
       ! reads lines until it reaches a defined material
       num_mat = 1
       do while (num_mat >= 0 .and. ier == 0)
          read(IIN_DB,'(A)',iostat=ier) line
          if (ier /= 0) exit

          ! skip empty/comment lines
          if (len_trim(line) == 0) cycle
          if (line(1:1) == '#' .or. line(1:1) == '!') cycle

          read(line,*) idomain_id,num_mat
       enddo
       if (ier /= 0) stop 'Error reading in undefined materials in nummaterial_velocity_file'

       ! checks if interface or tomography definition
       read(line,*) undef_mat_prop(6,imat),undef_mat_prop(1,imat),undef_mat_prop(2,imat)
       read(undef_mat_prop(1,imat),*) num_mat

       if (trim(undef_mat_prop(2,imat)) == 'interface') then
          ! line will have 5 arguments, e.g.: 2 -1 interface 1 2
          read(line,*) undef_mat_prop(6,imat),undef_mat_prop(1,imat),undef_mat_prop(2,imat), &
               undef_mat_prop(3,imat),undef_mat_prop(4,imat)
          undef_mat_prop(5,imat) = "0" ! dummy value
       else if (trim(undef_mat_prop(2,imat)) == 'tomography') then
          ! line will have 6 arguments, e.g.: 2 -1 tomography elastic tomography_model.xyz 1
          read(line,*) undef_mat_prop(6,imat),undef_mat_prop(1,imat),undef_mat_prop(2,imat), &
               undef_mat_prop(3,imat),undef_mat_prop(4,imat)
          undef_mat_prop(5,imat) = "0" ! dummy value
       else
          print *,'Error invalid type for undefined materials: ',trim(undef_mat_prop(2,imat))
          stop 'Error in nummaterial_velocity_file: invalid line in for undefined material'
       endif

       ! checks material_id
       if (trim(undef_mat_prop(2,imat)) == 'interface' .or. trim(undef_mat_prop(2,imat)) == 'tomography') then
          if (num_mat /= -imat) then
             print *,'Error undefined materials: invalid material id ',num_mat,' for undefined material',imat
             print *,'Please check line: ',trim(line)
             print *,'Note that material ids must be ordered in increasing (negative) order.'
             stop 'Error in nummaterial_velocity_file: Invalid material id for undefined materials.'
          endif
          if (num_mat > 0) &
               stop 'Error in nummaterial_velocity_file: Positive material id is invalid for undefined materials'
          if (-num_mat > count_undef_mat)  &
               stop 'Error in nummaterial_velocity_file: Negative material id exceeds bounds for undefined materials.'
       endif

       ! checks interface: flag_down/flag_up
       if (trim(undef_mat_prop(2,imat)) == 'interface') then
          ! flag_down
          read( undef_mat_prop(3,imat),*) num_mat
          if (num_mat > 0) then
             ! must point to a defined material
             if (num_mat > count_def_mat) &
                  stop "Error in nummaterial_velocity_file: invalid flag_down in interface definition"
          else
             ! must point to an undefined material
             if (-num_mat > count_undef_mat) &
                  stop "Error in nummaterial_velocity_file: invalid flag_down in interface definition"
          endif
          ! flag_up
          read( undef_mat_prop(4,imat),*) num_mat
          if (num_mat > 0) then
             ! must point to a defined material
             if (num_mat > count_def_mat) &
                  stop "Error in nummaterial_velocity_file: invalid flag_up in interface definition"
          else
             ! must point to an undefined material
             if (-num_mat > count_undef_mat) &
                  stop "Error in nummaterial_velocity_file: invalid flag_up in interface definition"
          endif
       endif

       ! checks domain id
       if (trim(undef_mat_prop(2,imat)) == 'tomography') then
          ! poroelastic domains not supported yet
          if (idomain_id /= 1 .and. idomain_id /= 2) &
          stop "Error in nummaterial_velocity_file: invalid domain id for undefined material, only acoustic/elastic supported."
          ! ignoring acoustic/elastic type name
          ! acoustic domains
          !if (idomain_id == 1 .and. trim(undef_mat_prop(3,imat)) /= 'acoustic') &
          !  stop "Error in nummaterial_velocity_file: invalid acoustic definition for undefined material"
          ! elastic domains
          !if (idomain_id == 2 .and. trim(undef_mat_prop(3,imat)) /= 'elastic') &
          !  stop "Error in nummaterial_velocity_file: invalid elastic definition for undefined material"
       endif
    enddo
    if (use_poroelastic_file) close(IIN_DB2)
    close(IIN_DB)

    do ispec=1,nspec_glob
       ! get material_id
       num_mat = mat(1,ispec)
       if (num_mat < 0) then
          ! finds undefined material property
          do imat=1,count_undef_mat
             if (-imat == num_mat) then
                ! interface
                if (trim(undef_mat_prop(2,imat)) == 'interface') then
                   mat(2,ispec) = 1
                   ! tomography
                else if (trim(undef_mat_prop(2,imat)) == 'tomography') then
                   mat(2,ispec) = 2
                else
                   ! shouldn't encounter this case
                   stop "Error in nummaterial_velocity_file: type name not recognized for undefined material"
                endif
             endif
          enddo
       else
          ! ispec belongs to a defined material
          mat(2,ispec) = 0
       endif
    enddo

    ! reads in absorbing boundary files
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_xmin', &
         status='old', form='formatted',iostat=ier)
    ! if the file does not exist then define the number of Stacey elements as zero for this face;
    ! beware that these files can also be used to set Dirichlet boundary conditions on the outer edges of CPML
    ! absorbing layers for elastic elements, not only for Stacey; thus these files may exist and be non-empty
    ! even when STACEY_ABSORBING_CONDITIONS is false
    if (ier /= 0) then
       nspec2D_xmin = 0
    else
       read(IIN_DB,*) nspec2D_xmin
    endif

    ! an array of size 0 is a valid object in Fortran 90, i.e. the array is then considered as allocated
    ! and can thus for instance be used as an argument in a call to a subroutine without giving any error
    ! even when full range and pointer checking is used in the compiler options;
    ! thus here the idea is that if some of the absorbing files do not exist because there are no absorbing
    ! conditions for this mesh then the array is created nonetheless, but with a dummy size of 0
    allocate(ibelm_xmin(nspec2D_xmin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 123')
    if (ier /= 0) stop 'Error allocating array ibelm_xmin'
    allocate(nodes_ibelm_xmin(NGNOD2D,nspec2D_xmin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 124')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_xmin'
    do ispec2D = 1,nspec2D_xmin
       ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
       ! note: ordering for CUBIT seems such that the normal of the face points outward of the element the face belongs to;
       !         in other words, nodes are in increasing order such that when looking from within the element outwards,
       !         they are ordered clockwise
       !
       !          doesn't necessarily have to start on top-rear, then bottom-rear, bottom-front, and finally top-front i.e.:
       !          point 1 = (0,1,1), point 2 = (0,1,0), point 3 = (0,0,0), point 4 = (0,0,1)
       read(IIN_DB,*) ibelm_xmin(ispec2D), (nodes_ibelm_xmin(inode,ispec2D), inode=1,NGNOD2D)
    enddo
    close(IIN_DB)
    write(27,*) 'absorbing boundaries:'
    write(27,*) '  nspec2D_xmin = ', nspec2D_xmin

    ! reads in absorbing boundary files
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_xmax', &
         status='old', form='formatted',iostat=ier)
    if (ier /= 0) then
       nspec2D_xmax = 0
    else
       read(IIN_DB,*) nspec2D_xmax
    endif
    allocate(ibelm_xmax(nspec2D_xmax),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 125')
    if (ier /= 0) stop 'Error allocating array ibelm_xmax'
    allocate(nodes_ibelm_xmax(NGNOD2D,nspec2D_xmax),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 126')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_xmax'
    do ispec2D = 1,nspec2D_xmax
       ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(IIN_DB,*) ibelm_xmax(ispec2D), (nodes_ibelm_xmax(inode,ispec2D), inode=1,NGNOD2D)
    enddo
    close(IIN_DB)
    write(27,*) '  nspec2D_xmax = ', nspec2D_xmax

  ! reads in absorbing boundary files
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_ymin', &
          status='old', form='formatted',iostat=ier)
    if (ier /= 0) then
       nspec2D_ymin = 0
    else
       read(IIN_DB,*) nspec2D_ymin
    endif
    allocate(ibelm_ymin(nspec2D_ymin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 127')
    if (ier /= 0) stop 'Error allocating array ibelm_ymin'
    allocate(nodes_ibelm_ymin(NGNOD2D,nspec2D_ymin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 128')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_ymin'
    do ispec2D = 1,nspec2D_ymin
       ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
       read(IIN_DB,*) ibelm_ymin(ispec2D), (nodes_ibelm_ymin(inode,ispec2D), inode=1,NGNOD2D)
    enddo
    close(IIN_DB)
    write(27,*) '  nspec2D_ymin = ', nspec2D_ymin

    ! reads in absorbing boundary files
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_ymax', &
         status='old', form='formatted',iostat=ier)
    if (ier /= 0) then
       nspec2D_ymax = 0
    else
       read(IIN_DB,*) nspec2D_ymax
    endif
    allocate(ibelm_ymax(nspec2D_ymax),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 129')
    if (ier /= 0) stop 'Error allocating array ibelm_ymax'
    allocate(nodes_ibelm_ymax(NGNOD2D,nspec2D_ymax),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 130')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_ymax'
    do ispec2D = 1,nspec2D_ymax
       ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
       read(IIN_DB,*) ibelm_ymax(ispec2D), (nodes_ibelm_ymax(inode,ispec2D), inode=1,NGNOD2D)
    enddo
    close(IIN_DB)
    write(27,*) '  nspec2D_ymax = ', nspec2D_ymax

    ! reads in absorbing boundary files
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_bottom', &
         status='old', form='formatted',iostat=ier)
    if (ier /= 0) then
       nspec2D_bottom = 0
    else
       read(IIN_DB,*) nspec2D_bottom
    endif
    allocate(ibelm_bottom(nspec2D_bottom),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 131')
    if (ier /= 0) stop 'Error allocating array ibelm_bottom'
    allocate(nodes_ibelm_bottom(NGNOD2D,nspec2D_bottom),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 132')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_bottom'
    do ispec2D = 1,nspec2D_bottom
       ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
       read(IIN_DB,*) ibelm_bottom(ispec2D), (nodes_ibelm_bottom(inode,ispec2D), inode=1,NGNOD2D)
    enddo
    close(IIN_DB)
    write(27,*) '  nspec2D_bottom = ', nspec2D_bottom

    ! reads in free_surface boundary files
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/free_or_absorbing_surface_file_zmax', &
         status='old', form='formatted',iostat=ier)
    if (ier /= 0) then
       nspec2D_top = 0
    else
       read(IIN_DB,*) nspec2D_top
    endif
    allocate(ibelm_top(nspec2D_top),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 133')
    if (ier /= 0) stop 'Error allocating array ibelm_top'
    allocate(nodes_ibelm_top(NGNOD2D,nspec2D_top),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 134')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_top'
    do ispec2D = 1,nspec2D_top
       ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
       read(IIN_DB,*) ibelm_top(ispec2D), (nodes_ibelm_top(inode,ispec2D), inode=1,NGNOD2D)
    enddo
    close(IIN_DB)
    write(27,*) '  nspec2D_top = ', nspec2D_top

    ! an array of size 0 is a valid object in Fortran 90, i.e. the array is then considered as allocated
    ! and can thus for instance be used as an argument in a call to a subroutine without giving any error
    ! even when full range and pointer checking is used in the compiler options;
    ! thus here the idea is that if some of the absorbing files do not exist because there are no absorbing
    ! conditions for this mesh then the array is created nonetheless, but with a dummy size of 0

    ! reads in absorbing_cpml boundary file
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_cpml_file', &
         status='old', form='formatted',iostat=ier)
    ! if the file does not exist but if there are PML_CONDITIONS then stop
    if (ier /= 0 .and. PML_CONDITIONS) &
         stop 'Error: PML_CONDITIONS is set to true but file absorbing_cpml_file does not exist'
    ! if the file does not exist or if there are no PML_CONDITIONS then define the number of CPML elements as zero
    if (ier /= 0 .or. .not. PML_CONDITIONS) then
       nspec_cpml = 0
    else
       read(IIN_DB,*) nspec_cpml
    endif

    ! sanity check
    if (PML_CONDITIONS .and. nspec_cpml <= 0) &
         stop 'Error: PML_CONDITIONS is set to true but nspec_cpml <= 0 in file absorbing_cpml_file'

    ! C-PML spectral elements global indexing
    allocate(cpml_to_spec(nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 135')
    if (ier /= 0) stop 'Error allocating array CPML_to_spec'
    ! C-PML regions (see below)
    allocate(cpml_regions(nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 136')
    if (ier /= 0) stop 'Error allocating array CPML_regions'
    do ispec_cpml=1,nspec_cpml
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
       read(IIN_DB,*) cpml_to_spec(ispec_cpml), cpml_regions(ispec_cpml)
    enddo
    close(IIN_DB)
    if (nspec_cpml > 0) write(27,*)  '  nspec_cpml = ', nspec_cpml

    ! sets mask of C-PML elements for all elements in this partition
    allocate(is_CPML(nspec_glob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 137')
    if (ier /= 0) stop 'Error allocating array is_CPML'
    is_CPML(:) = .false.
    do ispec_cpml=1,nspec_cpml
       if ((cpml_regions(ispec_cpml) >= 1) .and. (cpml_regions(ispec_cpml) <= 7)) then
          is_CPML(cpml_to_spec(ispec_cpml)) = .true.
       endif
    enddo

    ! reads in moho_surface boundary files (optional)
    open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/moho_surface_file', &
         status='old', form='formatted',iostat=ier)
    if (ier /= 0) then
       nspec2D_moho = 0
    else
       read(IIN_DB,*) nspec2D_moho
    endif
    allocate(ibelm_moho(nspec2D_moho),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 138')
    if (ier /= 0) stop 'Error allocating array ibelm_moho'
    allocate(nodes_ibelm_moho(NGNOD2D,nspec2D_moho),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 139')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_moho'
    do ispec2D = 1,nspec2D_moho
       ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
       read(IIN_DB,*) ibelm_moho(ispec2D), (nodes_ibelm_moho(inode,ispec2D), inode=1,NGNOD2D)
    enddo
    close(IIN_DB)
    if (nspec2D_moho > 0) write(27,*) '  nspec2D_moho = ', nspec2D_moho

    call read_fault_files(localpath_name)

    if (ANY_FAULT) then
       call save_nodes_coords(nodes_coords_glob,nnodes_glob)
       call close_faults(nodes_coords_glob,nnodes_glob)
    endif

  end subroutine read_mesh_files

  !--------------------------------------------------
  ! interface to loading : sets weights for acoustic/elastic/poroelastic elements to account for different
  !               expensive calculations in specfem simulations
  !--------------------------------------------------
  subroutine compute_load_elemnts()
    !!
    !! calling subroutine to fill load_elmnts array contains the element weight to
    !! be considered in mesh decomposition in order to obtain a good load inbalance
    !!
    implicit none

    integer :: ier

    allocate(load_elmnts(nspec_glob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 140')

    call  acoustic_elastic_poro_load (load_elmnts,nspec_glob,count_def_mat,count_undef_mat, &
                                    num_material,mat_prop,undef_mat_prop,ATTENUATION)

  end subroutine compute_load_elemnts


end module module_mesh


!-----------------------------  other subroutines -------------------------------------------------------

  !--------------------------------------------------
  ! loading : sets weights for acoustic/elastic/poroelastic elements to account for different
  !               expensive calculations in specfem simulations
  !--------------------------------------------------

subroutine acoustic_elastic_poro_load (elmnts_load,nspec,count_def_mat,count_undef_mat, &
                                    num_material,mat_prop,undef_mat_prop,ATTENUATION)

  use module_mesh, only: MAX_STRING_LEN, ACOUSTIC_LOAD, ELASTIC_LOAD, VISCOELASTIC_LOAD, POROELASTIC_LOAD
  !
  ! note:
  !   acoustic material    = domainID 1  (stored in mat_prop(6,..) )
  !   elastic material     = domainID 2
  !   poroelastic material = domainID 3
  !
  implicit none

  integer,                                                        intent(in)  :: nspec
  integer,                                                        intent(in)  :: count_def_mat,count_undef_mat
  logical,                                                        intent(in)  :: ATTENUATION

  ! materials
  integer,                       dimension(1:nspec),              intent(in)  :: num_material
  double precision,              dimension(16,count_def_mat),     intent(in)  :: mat_prop
  character(len=MAX_STRING_LEN), dimension(6,count_undef_mat),    intent(in)  :: undef_mat_prop

  ! load weights
  double precision,              dimension(1:nspec),              intent(out) :: elmnts_load

  ! local parameters
  logical, dimension(:), allocatable  :: is_acoustic, is_elastic, is_poroelastic
  integer  :: i,el,idomain_id,ier

  allocate(is_acoustic(-count_undef_mat:count_def_mat),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 141')
  allocate(is_elastic(-count_undef_mat:count_def_mat),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 142')
  allocate(is_poroelastic(-count_undef_mat:count_def_mat),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 143')

  ! initializes flags
  is_acoustic(:) = .false.
  is_elastic(:) = .false.
  is_poroelastic(:) = .false.


  ! sets acoustic/elastic/poroelastic flags for defined materials
  do i = 1, count_def_mat
     idomain_id = nint(mat_prop(7,i))
     ! acoustic material has idomain_id 1
     if (idomain_id == 1) then
        is_acoustic(i) = .true.
     endif
     ! elastic material has idomain_id 2
     if (idomain_id == 2) then
        is_elastic(i) = .true.
     endif
     ! poroelastic material has idomain_id 3
     if (idomain_id == 3) then
        is_poroelastic(i) = .true.
     endif
  enddo

  ! sets acoustic/elastic flags for undefined materials
  do i = 1, count_undef_mat
     read(undef_mat_prop(6,i),*) idomain_id
     ! acoustic material has idomain_id 1
     if (idomain_id == 1) then
        is_acoustic(-i) = .true.
     endif
     ! elastic material has idomain_id 2
     if (idomain_id == 2) then
        is_elastic(-i) = .true.
     endif
  enddo


  ! sets weights for elements
  do el = 0, nspec-1
     ! note: num_material index can be negative for tomographic material definitions
     ! acoustic element (cheap)
     if (num_material(el+1) > 0) then
        if (is_acoustic(num_material(el+1))) elmnts_load(el+1) = ACOUSTIC_LOAD
        ! elastic element (expensive)
        if (is_elastic(num_material(el+1))) then
           if (ATTENUATION) then
              elmnts_load(el+1) = VISCOELASTIC_LOAD
           else
              elmnts_load(el+1) = ELASTIC_LOAD
           endif
        endif
        ! poroelastic element (very expensive)
        if (is_poroelastic(num_material(el+1))) elmnts_load(el+1) = POROELASTIC_LOAD
     else ! JC JC: beware! To modify to take into account the -200? flags used in C-PML boundary conditions
        ! tomographic materials count as elastic
        if (ATTENUATION) then
           elmnts_load(el+1) = VISCOELASTIC_LOAD
        else
           elmnts_load(el+1) = ELASTIC_LOAD
        endif
     endif
  enddo

end subroutine acoustic_elastic_poro_load






