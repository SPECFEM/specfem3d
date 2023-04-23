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


  subroutine read_mesh_files()

! reads in mesh files
!
! note: this routine will be called by serial decompose mesher.

  use constants, only: NDIM,IIN_DB

  use decompose_mesh_par

  implicit none

  ! local parameters
  integer :: aniso_flag,idomain_id
  double precision :: vp,vs,rho,qkappa,qmu
  ! poroelastic parameters read in a new file
  double precision :: rhos,rhof,phi,tort,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,eta,mufr
  integer(kind=8) :: nspec_long
  integer :: inode,idummy
  integer :: ispec,ispec2D,ispec_CPML,ier

  character(len=MAX_STRING_LEN) :: line
  logical :: use_poroelastic_file
  logical :: file_found

  ! user output
  print *, 'reading mesh files in: ',trim(localpath_name)
  print *
  print *, '  using NGNOD = ',NGNOD
  if (NGNOD == 8) then
    print *, '  linear elements'
  else if (NGNOD == 27) then
    print *, '  quadratic elements'
  endif
  print *

  ! reads node coordinates
  open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/nodes_coords_file', &
        status='old', form='formatted', iostat = ier)
  if (ier /= 0) then
    print *,'could not open file:',localpath_name(1:len_trim(localpath_name))//'/nodes_coords_file'
    stop 'Error opening file nodes_coords_file'
  endif

  read(IIN_DB,*) nnodes

  if (nnodes < 1) stop 'Error: nnodes < 1'
  allocate(nodes_coords(NDIM,nnodes),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 85')
  if (ier /= 0) stop 'Error allocating array nodes_coords'
  nodes_coords(:,:) = 0.0d0
  do inode = 1, nnodes
    ! format: #id_node #x_coordinate #y_coordinate #z_coordinate
    read(IIN_DB,*) num_node, nodes_coords(1,num_node), nodes_coords(2,num_node), nodes_coords(3,num_node)

    ! bounds check
    if (num_node > nnodes .or. num_node < 1)  stop "Error : Invalid nodes_coords_file"
  enddo
  close(IIN_DB)
  print *, 'total number of nodes: '
  print *, '  nnodes = ', nnodes

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
  nspec = nspec_long

  if (nspec < 1) stop 'Error: nspec < 1'

  allocate(elmnts(NGNOD,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 86')
  if (ier /= 0) stop 'Error allocating array elmnts'
  elmnts(:,:) = 0
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
    !
    ! NGNOD == 8 ordering:
    !   SPECFEM ordering
    !           8 ------ 7
    !          /|      / |            z
    !         5 ----- 6  |            |
    !         | |     |  |            |
    !         | 4 ----|- 3            |-------> y
    !         |/      |/               \
    !         1 ----- 2                 \ x
    !
    ! NGNOD == 27 ordering:
    !   SPECFEM ordering
    !   1. corner nodes  1-8: same as hex8
    !   2. midside nodes (nodes located in the middle of an edge): edges bottom, mid (bottom-to-top), top
    !   3. side center nodes (nodes located in the middle of a face): bottom,front,right,back,left,top
    !   4. center node
    !
    !      8----19-----7
    !      |\          |\
    !      | 20    26  | 18
    !      16  \ 24    15 \
    !      |    5----17+---6
    !      | 25 |  27  | 23|
    !      4----+-11---3   |
    !       \  13    22 \  14
    !        12 |  21    10|
    !          \|         \|
    !           1----9----2
    !
    ! compare with http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    !      3----13-----2
    !      |\          |\
    !      | 15    24  | 14        v
    !      9  \ 20     11 \        |
    !      |    7----19+---6       |
    !      | 22 |  26  | 23|       |-----> u
    !      0----+-8----1   |        \
    !       \  17    25 \  18        \
    !        10 |  21    12|          w
    !          \|         \|
    !           4----16----5
    !
    read(IIN_DB,*,iostat=ier) num_elmnt,(elmnts(inode,num_elmnt), inode=1,NGNOD)

    ! debug
    !print *,'ispec: ',ispec,'num_elmnt:',num_elmnt,elmnts(1:NGNOD,num_elmnt)

    if (ier /= 0) then
      print *,'Error while attempting to read ',NGNOD,'element data values from the mesh file'
      print *,'  ispec: ',ispec,'num_elmnt: ',num_elmnt,'elemnts',elmnts(1:NGNOD,num_elmnt)
      if (NGNOD == 8) print *,'check if your mesh file is indeed composed of HEX8 elements'
      if (NGNOD == 27) print *,'check if your mesh file is indeed composed of HEX27 elements'
      stop 'Error reading element data from the mesh file'
    endif

    ! bounds check
    if (num_elmnt > nspec .or. num_elmnt < 1)  stop "Error : Invalid mesh_file"

  enddo
  close(IIN_DB)
  print *, 'total number of spectral elements:'
  print *, '  nspec = ', nspec

  ! reads material associations
  open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/materials_file', &
        status='old', form='formatted',iostat=ier)
  if (ier /= 0) stop 'Error opening materials_file'

  allocate(mat(2,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 87')
  if (ier /= 0) stop 'Error allocating array mat'
  mat(:,:) = 0

  do ispec = 1, nspec
    ! format: #id_element #flag
    ! note: be aware that elements may not be sorted in materials_file
    read(IIN_DB,*) num_mat,mat(1,num_mat)
    if ((num_mat > nspec) .or. (num_mat < 1)) stop "Error : Invalid materials_file"
  enddo
  close(IIN_DB)

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
  print *,'materials:'
  ier = 0
  do while (ier == 0)
    ! note: format #material_domain_id #material_id #...
    read(IIN_DB,'(A)',iostat=ier) line
    if (ier /= 0) exit

    ! skip empty/comment lines
    line = adjustl(line)        ! suppress leading white spaces, if any
    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#' .or. line(1:1) == '!') cycle

    read(line,*,iostat=ier) idummy,num_mat
    if (ier /= 0) exit

    print *, '  num_mat = ',num_mat
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
  print *, '  defined = ',count_def_mat, 'undefined = ',count_undef_mat

  ! check with material flags
  ! defined materials
  if (count_def_mat > 0 .and. maxval(mat(1,:)) > count_def_mat) then
    print *,'Error material definitions:'
    print *,'  materials associated in materials_file:',maxval(mat(1,:))
    print *,'  larger than defined materials in nummaterial_velocity_file:',count_def_mat
    stop 'Error positive material id exceeds bounds for defined materials'
  endif

  ! material properties
  allocate(mat_prop(17,count_def_mat),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 88')
  if (ier /= 0) stop 'Error allocating array mat_prop'
  mat_prop(:,:) = 0.d0

  ! undefined materials
  if (count_undef_mat > 0 .and. minval(mat(1,:)) < -count_undef_mat) then
    print *,'Error material definitions:'
    print *,'  materials associated in materials_file:',minval(mat(1,:))
    print *,'  larger than undefined materials in nummaterial_velocity_file:',count_undef_mat
    stop 'Error negative material id exceeds bounds for undefined materials'
  endif

  ! undefined material properties (will get determined from tomographic model properties)
  allocate(undef_mat_prop(6,count_undef_mat),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 89')
  if (ier /= 0) stop 'Error allocating array undef_mat_prop'
  undef_mat_prop(:,:) = ''

  ! reads in defined material properties
  rewind(IIN_DB,iostat=ier)
  if (ier /= 0) stop 'Error rewinding nummaterial_velocity_file'

  ! poroelastic materials
  !
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
  if (ier /= 0) then
    use_poroelastic_file = .false.
    !stop 'Error opening nummaterial_poroelastic_file'
    print *, '  no poroelastic material file found'
  else
    use_poroelastic_file = .true.
    print *, '  poroelastic material file found'
  endif
  ier = 0

  ! note: entries in nummaterial_velocity_file can be an unsorted list of all
  !          defined materials (material_id > 0) and undefined materials (material_id < 0)
  print *, '  defined materials  : ',count_def_mat
  do imat = 1,count_def_mat
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
      line = adjustl(line)        ! suppress leading white spaces, if any
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle

      read(line,*) idomain_id,num_mat
    enddo
    if (ier /= 0) stop 'Error reading in defined materials in nummaterial_velocity_file'

    ! reads in defined material properties
    read(line,*,iostat=ier) idomain_id,num_mat,rho,vp,vs,qkappa,qmu,aniso_flag
    if (ier /= 0) stop 'Error reading nummaterial_velocity_file, please check if it has the required format...'

    ! sanity check: Q factor cannot be equal to zero, thus convert to 9999 to indicate no attenuation
    ! if users have used 0 to indicate that instead
    if (qkappa <= 0.000001) qkappa = 9999.d0
    if (qmu <= 0.000001) qmu = 9999.d0

    ! checks material_id bounds
    if (num_mat < 1 .or. num_mat > count_def_mat) &
      stop "Error in nummaterial_velocity_file: material id invalid for defined materials."

    ! reads in material properties
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
      do while (ier == 0)
        read(97,'(A)',iostat=ier) line
        if (ier /= 0) exit

        ! skip empty/comment lines
        line = adjustl(line)        ! suppress leading white spaces, if any
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#' .or. line(1:1) == '!') cycle

        ! format:
        ! #rho_s #rho_f #phi #tort #kxx #kxy #kxz #kyy #kyz #kzz #kappa_s #kappa_f #kappa_fr #eta #mu_fr
        read(line,*,iostat=ier) rhos,rhof,phi,tort,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,eta,mufr
        if (ier /= 0) then
          stop 'Error reading nummaterial_poroelastic_file, please check if it has the required format...'
        endif

        ! only read 1 valid line at a time
        if (ier == 0) exit
      enddo
      if (ier /= 0) stop 'Error reading in materials in nummaterial_poroelastic_file'

      mat_prop(1,num_mat) = rhos
      mat_prop(2,num_mat) = rhof
      mat_prop(3,num_mat) = phi
      mat_prop(4,num_mat) = tort
      mat_prop(5,num_mat) = eta
      mat_prop(6,num_mat) = 0       ! aniso_flag
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
  print *, '  undefined materials: ',count_def_mat,' (interfaces/tomography models/..)'
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
       line = adjustl(line)        ! suppress leading white spaces, if any
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#' .or. line(1:1) == '!') cycle

       read(line,*,iostat=ier) idomain_id,num_mat
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
  if (use_poroelastic_file) close(97)
  close(IIN_DB)

  do ispec = 1,nspec
    ! get material_id
    num_mat = mat(1,ispec)
    if (num_mat < 0) then
      ! finds undefined material property
      do imat = 1,count_undef_mat
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
    print *, '  no absorbing_surface_file_xmin file found'
  else
    read(IIN_DB,*) nspec2D_xmin
  endif

! an array of size 0 is a valid object in Fortran 90, i.e. the array is then considered as allocated
! and can thus for instance be used as an argument in a call to a subroutine without giving any error
! even when full range and pointer checking is used in the compiler options;
! thus here the idea is that if some of the absorbing files do not exist because there are no absorbing
! conditions for this mesh then the array is created nonetheless, but with a dummy size of 0

  if (nspec2D_xmin > 0) then
    allocate(ibelm_xmin(nspec2D_xmin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 90')
    if (ier /= 0) stop 'Error allocating array ibelm_xmin'
    allocate(nodes_ibelm_xmin(NGNOD2D,nspec2D_xmin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 91')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_xmin'
  else
    ! dummy allocation
    allocate(ibelm_xmin(1),nodes_ibelm_xmin(1,1))
  endif
  ibelm_xmin(:) = 0; nodes_ibelm_xmin(:,:) = 0

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
  print *, 'absorbing boundaries:'
  print *, '  nspec2D_xmin = ', nspec2D_xmin

  ! reads in absorbing boundary files
  open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_xmax', &
        status='old', form='formatted',iostat=ier)
  if (ier /= 0) then
    nspec2D_xmax = 0
    print *, '  no absorbing_surface_file_xmax file found'
  else
    read(IIN_DB,*) nspec2D_xmax
  endif

  if (nspec2D_xmax > 0) then
    allocate(ibelm_xmax(nspec2D_xmax),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 92')
    if (ier /= 0) stop 'Error allocating array ibelm_xmax'
    allocate(nodes_ibelm_xmax(NGNOD2D,nspec2D_xmax),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 93')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_xmax'
  else
    ! dummy allocation
    allocate(ibelm_xmax(1),nodes_ibelm_xmax(1,1))
  endif
  ibelm_xmax(:) = 0; nodes_ibelm_xmax(:,:) = 0

  do ispec2D = 1,nspec2D_xmax
    ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
    read(IIN_DB,*) ibelm_xmax(ispec2D), (nodes_ibelm_xmax(inode,ispec2D), inode=1,NGNOD2D)
  enddo
  close(IIN_DB)
  print *, '  nspec2D_xmax = ', nspec2D_xmax

  ! reads in absorbing boundary files
  open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_ymin', &
        status='old', form='formatted',iostat=ier)
  if (ier /= 0) then
    nspec2D_ymin = 0
    print *, '  no absorbing_surface_file_ymin file found'
  else
    read(IIN_DB,*) nspec2D_ymin
  endif

  if (nspec2D_ymin > 0) then
    allocate(ibelm_ymin(nspec2D_ymin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 94')
    if (ier /= 0) stop 'Error allocating array ibelm_ymin'
    allocate(nodes_ibelm_ymin(NGNOD2D,nspec2D_ymin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 95')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_ymin'
  else
    ! dummy allocation
    allocate(ibelm_ymin(1),nodes_ibelm_ymin(1,1))
  endif
  ibelm_ymin(:) = 0; nodes_ibelm_ymin(:,:) = 0

  do ispec2D = 1,nspec2D_ymin
    ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
    read(IIN_DB,*) ibelm_ymin(ispec2D), (nodes_ibelm_ymin(inode,ispec2D), inode=1,NGNOD2D)
  enddo
  close(IIN_DB)
  print *, '  nspec2D_ymin = ', nspec2D_ymin

  ! reads in absorbing boundary files
  open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_ymax', &
        status='old', form='formatted',iostat=ier)
  if (ier /= 0) then
    nspec2D_ymax = 0
    print *, '  no absorbing_surface_file_ymax file found'
  else
    read(IIN_DB,*) nspec2D_ymax
  endif

  if (nspec2D_ymax > 0) then
    allocate(ibelm_ymax(nspec2D_ymax),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 96')
    if (ier /= 0) stop 'Error allocating array ibelm_ymax'
    allocate(nodes_ibelm_ymax(NGNOD2D,nspec2D_ymax),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 97')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_ymax'
  else
    ! dummy allocation
    allocate(ibelm_ymax(1),nodes_ibelm_ymax(1,1))
  endif
  ibelm_ymax(:) = 0; nodes_ibelm_ymax(:,:) = 0

  do ispec2D = 1,nspec2D_ymax
    ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
    read(IIN_DB,*) ibelm_ymax(ispec2D), (nodes_ibelm_ymax(inode,ispec2D), inode=1,NGNOD2D)
  enddo
  close(IIN_DB)
  print *, '  nspec2D_ymax = ', nspec2D_ymax

  ! reads in absorbing boundary files
  open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_bottom', &
        status='old', form='formatted',iostat=ier)
  if (ier /= 0) then
    nspec2D_bottom = 0
    print *, '  no absorbing_surface_file_bottom file found'
  else
    read(IIN_DB,*) nspec2D_bottom
  endif

  if (nspec2D_bottom > 0) then
    allocate(ibelm_bottom(nspec2D_bottom),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 98')
    if (ier /= 0) stop 'Error allocating array ibelm_bottom'
    allocate(nodes_ibelm_bottom(NGNOD2D,nspec2D_bottom),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 99')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_bottom'
  else
    ! dummy allocation
    allocate(ibelm_bottom(1),nodes_ibelm_bottom(1,1))
  endif
  ibelm_bottom(:) = 0; nodes_ibelm_bottom(:,:) = 0

  do ispec2D = 1,nspec2D_bottom
    ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
    read(IIN_DB,*) ibelm_bottom(ispec2D), (nodes_ibelm_bottom(inode,ispec2D), inode=1,NGNOD2D)
  enddo
  close(IIN_DB)
  print *, '  nspec2D_bottom = ', nspec2D_bottom

  ! reads in free_surface boundary files
  ! checks old (now obsolete) naming formats to avoid mistakenly using obsolete files
  open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/free_surface_file', &
      status='old', form='formatted',iostat=ier)
  if (ier == 0) then
    print *, '  "free_surface_file" file name is deprecated!'
    print *, '  Please rename it to "free_or_absorbing_surface_file_zmax" and restart the run'
    print *, '  or remove it if it is obsolete'
    stop 'please rename or remove the input file with obsolete file name'
  endif

  open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_zmax', &
      status='old', form='formatted',iostat=ier)
  if (ier == 0) then
    print *, '  "absorbing_surface_file_zmax" file name is deprecated!'
    print *, '  Please rename it to "free_or_absorbing_surface_file_zmax" and restart the run'
    print *, '  or remove it if it is obsolete'
    stop 'please rename or remove the input file with obsolete file name'
  endif

  open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_top', &
      status='old', form='formatted',iostat=ier)
  if (ier == 0) then
    print *, '  "absorbing_surface_file_top" file name is deprecated!'
    print *, '  Please rename it to "free_or_absorbing_surface_file_zmax" and restart the run'
    print *, '  or remove it if it is obsolete'
    stop 'please rename or remove the input file with obsolete file name'
  endif

  file_found = .false.
  open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/free_or_absorbing_surface_file_zmax', &
        status='old', form='formatted',iostat=ier)
  if (ier == 0) file_found = .true.
  if (.not. file_found) then
    nspec2D_top = 0
    print *, '  no free_or_absorbing_surface_file_zmax file found'
  else
    read(IIN_DB,*) nspec2D_top
  endif

  if (nspec2D_top > 0) then
    allocate(ibelm_top(nspec2D_top),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 100')
    if (ier /= 0) stop 'Error allocating array ibelm_top'
    allocate(nodes_ibelm_top(NGNOD2D,nspec2D_top),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 101')
    if (ier /= 0) stop 'Error allocating array nodes_ibelm_top'
  else
    ! dummy allocation
    allocate(ibelm_top(1),nodes_ibelm_top(1,1))
  endif
  ibelm_top(:) = 0; nodes_ibelm_top(:,:) = 0

  do ispec2D = 1,nspec2D_top
    ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
    read(IIN_DB,*) ibelm_top(ispec2D), (nodes_ibelm_top(inode,ispec2D), inode=1,NGNOD2D)
  enddo
  close(IIN_DB)
  print *, '  nspec2D_top = ', nspec2D_top

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
  ! note: in case there is a cpml file, we will read it in and let the user decide when running an actual simulation
  !       if he wants to use cpml elements or not... otherwise he will need to re-run the whole meshing part when
  !       switching PML_CONDITIONS
  if (ier /= 0) then
    nspec_cpml = 0
    print *, '  no absorbing_cpml_file file found'
  else
    read(IIN_DB,*) nspec_cpml
  endif

  ! sanity check
  if (PML_CONDITIONS .and. nspec_cpml <= 0) &
      stop 'Error: PML_CONDITIONS is set to true but nspec_cpml <= 0 in file absorbing_cpml_file'

  ! C-PML spectral elements global indexing
  if (nspec_cpml > 0) then
    allocate(CPML_to_spec(nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 102')
    if (ier /= 0) stop 'Error allocating array CPML_to_spec'
    ! C-PML regions (see below)
    allocate(CPML_regions(nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 103')
    if (ier /= 0) stop 'Error allocating array CPML_regions'
  else
    ! dummy allocation
    allocate(CPML_to_spec(1),CPML_regions(1))
  endif
  CPML_to_spec(:) = 0; CPML_regions(:) = 0

  do ispec_CPML = 1,nspec_cpml
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
     read(IIN_DB,*) CPML_to_spec(ispec_CPML), CPML_regions(ispec_CPML)
  enddo
  close(IIN_DB)
  if (nspec_cpml > 0) print *, '  nspec_cpml = ', nspec_cpml

  ! sets mask of C-PML elements for all elements in this partition
  allocate(is_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 104')
  if (ier /= 0) stop 'Error allocating array is_CPML'
  is_CPML(:) = .false.

  do ispec_CPML = 1,nspec_cpml
     if ((CPML_regions(ispec_CPML) >= 1) .and. (CPML_regions(ispec_CPML) <= 7)) then
        is_CPML(CPML_to_spec(ispec_CPML)) = .true.
     endif
  enddo

  ! reads in moho_surface boundary files (optional)
  open(unit=IIN_DB, file=localpath_name(1:len_trim(localpath_name))//'/moho_surface_file', &
        status='old', form='formatted',iostat=ier)
  if (ier /= 0) then
    nspec2D_moho = 0
    print *, '  no moho_surface_file file found'
  else
    read(IIN_DB,*) nspec2D_moho
  endif
  allocate(ibelm_moho(nspec2D_moho),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 105')
  if (ier /= 0) stop 'Error allocating array ibelm_moho'
  allocate(nodes_ibelm_moho(NGNOD2D,nspec2D_moho),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 106')
  if (ier /= 0) stop 'Error allocating array nodes_ibelm_moho'
  do ispec2D = 1,nspec2D_moho
    ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
    read(IIN_DB,*) ibelm_moho(ispec2D), (nodes_ibelm_moho(inode,ispec2D), inode=1,NGNOD2D)
  enddo
  close(IIN_DB)
  if (nspec2D_moho > 0) print *, '  nspec2D_moho = ', nspec2D_moho

  ! fault surfaces
  call read_fault_files(localpath_name)
  if (ANY_FAULT) then
    ! saving original node coordinates, where split nodes are still open
    call save_nodes_coords(nodes_coords,nnodes)
    ! closing split node gaps in nodes_coords to have a conforming mesh
    call close_faults(nodes_coords,nnodes)
  endif

  end subroutine read_mesh_files

