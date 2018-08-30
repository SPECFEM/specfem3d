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



module create_meshfem_par

  use constants, only: MAX_STRING_LEN

  implicit none

  ! node coordinates
  double precision, dimension(:,:), allocatable :: nodes_coords

  ! material id
  integer, dimension(:), allocatable :: ispec_material_id

  ! boundary locator
  logical, dimension(:,:), allocatable :: iboun

  ! MPI cut-planes parameters along xi and along eta
  logical, dimension(:,:), allocatable :: iMPIcut_xi,iMPIcut_eta

  ! name of the database file
  character(len=MAX_STRING_LEN) :: prname

end module create_meshfem_par

!
!----------------------------------------------------------------------------------------------------
!


  subroutine create_meshfem_mesh()

  use create_meshfem_par, only: prname,nodes_coords,ispec_material_id,iboun,iMPIcut_xi,iMPIcut_eta

  use meshfem3D_par, only: NSPEC_AB,xgrid,ygrid,zgrid,ibool, &
    xstore,ystore,zstore, &
    iproc_xi_current,iproc_eta_current,addressing,nspec, &
    NGLOB_AB,npointot, &
    NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER, &
    NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    NPROC_XI,NPROC_ETA, &
    nsubregions,subregions,NMATERIALS,material_properties, &
    myrank, sizeprocs, &
    LOCAL_PATH,UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK,NEX_XI,NEX_ETA, &
    CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES, &
    USE_REGULAR_MESH,NDOUBLINGS,ner_doublings, &
    ADIOS_ENABLED, ADIOS_FOR_DATABASES, &
    nspec_CPML,is_CPML,CPML_to_spec,CPML_regions

  ! create the different regions of the mesh
  use constants, only: MF_IN_DATA_FILES,MAX_STRING_LEN,NGNOD_EIGHT_CORNERS, &
    IMAIN,IOVTK,CUSTOM_REAL,HUGEVAL,TINYVAL

  use constants_meshfem3D, only: NSPEC_DOUBLING_SUPERBRICK,NGLOB_DOUBLING_SUPERBRICK, &
    IFLAG_BASEMENT_TOPO,IFLAG_ONE_LAYER_TOPOGRAPHY, &
    NGLLCUBE_M,NGLLX_M,NGLLY_M,NGLLZ_M

  use adios_manager_mod, only: adios_setup,adios_cleanup

  implicit none

  ! local parameters
  ! auxiliary variables to generate the mesh
  integer :: isubregion,doubling_index,nmeshregions
  integer, dimension(:,:,:), allocatable :: material_num

  ! topology of the elements
  integer :: iaddx(NGNOD_EIGHT_CORNERS)
  integer :: iaddy(NGNOD_EIGHT_CORNERS)
  integer :: iaddz(NGNOD_EIGHT_CORNERS)

  double precision :: xelm(NGNOD_EIGHT_CORNERS)
  double precision :: yelm(NGNOD_EIGHT_CORNERS)
  double precision :: zelm(NGNOD_EIGHT_CORNERS)

  ! variables for creating array ibool (some arrays also used for AVS or DX files)
  integer, dimension(:), allocatable :: iglob,locval
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,zp

  integer :: nglob

  ! boundary parameters locator
  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top

  ! number of elements on the boundaries
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  ! flag indicating whether point is in the sediments
  !  logical point_is_in_sediments
  logical, dimension(:,:,:,:), allocatable :: flag_sediments
  logical, dimension(:), allocatable :: not_fully_in_bedrock

  ! material properties
  double precision :: VP_MAX

  ! doublings zone
  integer :: nspec_sb
  logical, dimension(NSPEC_DOUBLING_SUPERBRICK,6) :: iboun_sb
  integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick
  double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick

  ! **************

  !--- Initialize ADIOS and setup the buffer size
  if (ADIOS_ENABLED) then
    call adios_setup()
  endif

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

  ! allocates arrays
  call cmm_allocate_arrays()

  ! generate mesh elements
  call cmm_create_mesh_elements()

  ! creates indirect addressing
  call cmm_create_addressing()

  ! determine cavity
  call cmm_determine_cavity(nglob)

  ! CPML initialization
  call create_CPML_regions(nspec,nglob,nodes_coords)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'saving mesh files'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! outputs mesh file for visualization
  call create_visual_files(CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES, &
                          nspec,nglob, &
                          prname,nodes_coords,ibool,ispec_material_id)

  ! stores boundary informations
  call store_boundaries(myrank,iboun,nspec, &
                       ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                       nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                       NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

  ! checks mesh resolution
  VP_MAX = maxval(material_properties(:,2))
  call check_mesh_quality(myrank,VP_MAX,nglob,nspec, &
                        nodes_coords(:,1),nodes_coords(:,2),nodes_coords(:,3),ibool, &
                        CREATE_VTK_FILES,prname)


  ! saves mesh as databases file
  if (ADIOS_FOR_DATABASES) then
    call save_databases_adios(LOCAL_PATH, myrank, sizeprocs, &
                              nspec,nglob,iproc_xi_current,iproc_eta_current, &
                              NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta, &
                              ibool,nodes_coords,ispec_material_id, &
                              nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                              NSPEC2D_BOTTOM,NSPEC2D_TOP, NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                              ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                              NMATERIALS,material_properties, &
                              nspec_CPML,CPML_to_spec,CPML_regions,is_CPML)
  else
    ! saves mesh as databases file  !! VM VM added xstore, ystore, zstore used for Axisem Coupling
    call save_databases(prname,nspec,nglob,iproc_xi_current,iproc_eta_current, &
                        NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta, &
                        ibool,nodes_coords,ispec_material_id, &
                        nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                        NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                        ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                        NMATERIALS,material_properties, &
                        nspec_cpml,CPML_to_spec,CPML_regions,is_CPML, &
                        xstore, ystore, zstore)
  endif

  !--- Clean ADIOS. Make sure everything is already written
  if (ADIOS_ENABLED) then
    call adios_cleanup()
  endif

  ! frees memory
  deallocate(ispec_material_id)
  deallocate(locval,ifseg,xp,yp,zp)
  deallocate(nodes_coords)
  deallocate(ibool,xstore,ystore,zstore)
  deallocate(is_CPML)
  deallocate(CPML_to_spec)
  deallocate(CPML_regions)

contains

!
!----------------------------------------------------------------------------------------------------
!

    subroutine cmm_allocate_arrays()

    implicit none
    ! local parameters
    integer :: ier

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'allocating mesh arrays'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! assign theoretical number of elements
    nspec = NSPEC_AB

    ! compute maximum number of points
    npointot = nspec * NGLLCUBE_M

    ! make sure everybody is synchronized
    call synchronize_all()

    ! use dynamic allocation to allocate memory for arrays
    allocate(ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1257')
    if (ier /= 0) stop 'Error allocating array ibool'

    allocate(xstore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1258')
    if (ier /= 0) stop 'Error allocating array xstore'
    allocate(ystore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1259')
    if (ier /= 0) stop 'Error allocating array ystore'
    allocate(zstore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1260')
    ! exit if there is not enough memory to allocate all the arrays
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

    ! flag indicating whether point is in the sediments
    allocate(flag_sediments(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1261')
    if (ier /= 0) stop 'Error allocating array flag_sediments'
    allocate(not_fully_in_bedrock(nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1262')
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

    ! boundary locator
    allocate(iboun(6,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1263')
    if (ier /= 0) stop 'Error allocating array iboun'

    ! boundary parameters locator
    allocate(ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1264')
    if (ier /= 0) stop 'Error allocating array ibelm_xmin'
    allocate(ibelm_xmax(NSPEC2DMAX_XMIN_XMAX),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1265')
    if (ier /= 0) stop 'Error allocating array ibelm_xmax'
    allocate(ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1266')
    if (ier /= 0) stop 'Error allocating array ibelm_ymin'
    allocate(ibelm_ymax(NSPEC2DMAX_YMIN_YMAX),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1267')
    if (ier /= 0) stop 'Error allocating array ibelm_ymax'
    allocate(ibelm_bottom(NSPEC2D_BOTTOM),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1268')
    if (ier /= 0) stop 'Error allocating array ibelm_bottom'
    allocate(ibelm_top(NSPEC2D_TOP),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1269')
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

    ! MPI cut-planes parameters along xi and along eta
    allocate(iMPIcut_xi(2,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1270')
    if (ier /= 0) stop 'Error allocating array iMPIcut_xi'
    allocate(iMPIcut_eta(2,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1271')
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

    ! allocate memory for arrays
    allocate(iglob(npointot),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1272')
    if (ier /= 0) stop 'Error allocating array iglob'
    allocate(locval(npointot),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1273')
    if (ier /= 0) stop 'Error allocating array locval'
    allocate(ifseg(npointot),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1274')
    if (ier /= 0) stop 'Error allocating array ifseg'
    allocate(xp(npointot),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1275')
    if (ier /= 0) stop 'Error allocating array xp'
    allocate(yp(npointot),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1276')
    if (ier /= 0) stop 'Error allocating array yp'
    allocate(zp(npointot),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1277')
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

    ! allocate material ids array
    allocate(material_num(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1278')
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

    allocate(ispec_material_id(nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1279')
    if (ier /= 0) stop 'Error allocating array ispec_material_id'

    ! synchronize
    call synchronize_all()

    end subroutine cmm_allocate_arrays

!
!----------------------------------------------------------------------------------------------------
!

    subroutine cmm_create_mesh_elements()

    implicit none
    ! local parameters
    integer :: ix,iy,ir,ir1,ir2,dir
    integer :: ix1,ix2,dix,iy1,iy2,diy
    integer :: iax,iay,iar,ia
    integer :: imaterial_number
    integer :: ioffset_x,ioffset_y,ioffset_z
    integer :: ispec,ispec_superbrick
    double precision :: x_min,x_max,y_min,y_max,z_min,z_max
    double precision :: x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob

    ! generate the elements in all the regions of the mesh
    ispec = 0

    ! to check that we assign a material to each element
    material_num(:,:,:) = -1000 ! dummy value
    ispec_material_id(:) = -1000

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'number of subregions = ',nsubregions
      call flush_IMAIN()
    endif

    do isubregion = 1,nsubregions
      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  defining subregion ',isubregion
        call flush_IMAIN()
      endif

      call define_model_regions(NEX_PER_PROC_XI,NEX_PER_PROC_ETA,iproc_xi_current,iproc_eta_current, &
                                isubregion,nsubregions,subregions, &
                                iaddx,iaddy,iaddz,ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
                                imaterial_number)

      ! loop on all the mesh points in current subregion
      do ir = ir1,ir2,dir
        do iy = iy1,iy2,diy
          do ix = ix1,ix2,dix
            material_num(ir,ix,iy) = imaterial_number
          enddo
        enddo
      enddo

    enddo

    if (USE_REGULAR_MESH) then
      nmeshregions = 1
    else
      call define_superbrick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb)
      nspec_sb = NSPEC_DOUBLING_SUPERBRICK

      ! mesh subregion:
      !   example NDOUBLINGS = 1: isubregion = 1,2,3
      !   example NDOUBLINGS = 2: isubregion = 1,2,3,4,5
      ! the values goes from 1 to 2 * NDOUBLINGS + 1
      nmeshregions = 2 * NDOUBLINGS + 1

      ! hard-coded way
      !if (NDOUBLINGS == 1) then
      !  nmeshregions = 3
      !else if (NDOUBLINGS == 2) then
      !  nmeshregions = 5
      !else
      !  stop 'Wrong number of mesh doubling regions'
      !endif
    endif

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'number of mesh regions = ',nmeshregions
      call flush_IMAIN()
    endif

    do isubregion = 1,nmeshregions
      ! user output
      if (myrank == 0) then
        if (modulo(isubregion,2) == 1) then
          write(IMAIN,*) '  creating mesh region ',isubregion
        else
          write(IMAIN,*) '  creating mesh region ',isubregion,' with doubling layer'
        endif
        call flush_IMAIN()
      endif

      ! define shape of elements
      call define_mesh_regions(myrank,USE_REGULAR_MESH,isubregion,NER,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                               iproc_xi_current,iproc_eta_current, &
                               NDOUBLINGS,ner_doublings, &
                               iaddx,iaddy,iaddz,ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar)

      ! loop on all the mesh points in current subregion
      if (modulo(isubregion,2) == 1) then

        ! Regular subregion case

        do ir = ir1,ir2,dir
          do iy = iy1,iy2,diy
            do ix = ix1,ix2,dix
              ! loop over the NGNOD_EIGHT_CORNERS nodes
              do ia = 1,NGNOD_EIGHT_CORNERS
                ! define topological coordinates of this mesh point
                ioffset_x = ix+iax*iaddx(ia)
                ioffset_y = iy+iay*iaddy(ia)
                ioffset_z = ir+iar*iaddz(ia)

                xelm(ia) = xgrid(ioffset_z,ioffset_x,ioffset_y)
                yelm(ia) = ygrid(ioffset_z,ioffset_x,ioffset_y)
                zelm(ia) = zgrid(ioffset_z,ioffset_x,ioffset_y)
              enddo

              ! add one spectral element to the list and store its material number
              ispec = ispec + 1

              if (ispec > nspec) then
                call exit_MPI(myrank,'ispec greater than nspec in mesh creation')
              endif

              ! check if element is on topography
              if ((ir == ir2) .and. (isubregion == nmeshregions)) then
                ! sets index to topography layer
                doubling_index = IFLAG_ONE_LAYER_TOPOGRAPHY
              else
                ! sets index to deeper, basement layers
                doubling_index = IFLAG_BASEMENT_TOPO
              endif

              ispec_material_id(ispec) = material_num(ir,ix,iy)

              ! store coordinates
              call store_coords(xstore,ystore,zstore,xelm,yelm,zelm,ispec,nspec)

              ! detect mesh boundaries
              call get_flags_boundaries(nspec,iproc_xi_current,iproc_eta_current, &
                                        ispec,doubling_index, &
                                        xstore(:,:,:,ispec),ystore(:,:,:,ispec),zstore(:,:,:,ispec), &
                                        iboun,iMPIcut_xi,iMPIcut_eta,NPROC_XI,NPROC_ETA, &
                                        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK,NEX_XI,NEX_ETA)
            enddo
          enddo
        enddo
        ! end of loop on all the mesh points in current subregion

      else

        ! Irregular subregion case

        do ir = ir1,ir2,dir
          do iy = iy1,iy2,diy
            do ix = ix1,ix2,dix
              ! loop on all the elements in the mesh doubling superbrick
              do ispec_superbrick = 1,nspec_sb
                ! loop on all the corner nodes of this element
                do ia = 1,NGNOD_EIGHT_CORNERS
                  ! define topological coordinates of this mesh point
                  ioffset_x = int(ix + iax*x_superbrick(ibool_superbrick(ia,ispec_superbrick)))
                  ioffset_y = int(iy + iay*y_superbrick(ibool_superbrick(ia,ispec_superbrick)))
                  ioffset_z = int(ir + iar*z_superbrick(ibool_superbrick(ia,ispec_superbrick)))

                  xelm(ia) = xgrid(ioffset_z,ioffset_x,ioffset_y)
                  yelm(ia) = ygrid(ioffset_z,ioffset_x,ioffset_y)
                  zelm(ia) = zgrid(ioffset_z,ioffset_x,ioffset_y)
                enddo

                ! add one spectral element to the list and store its material number
                ispec = ispec + 1

                if (ispec > nspec) then
                  call exit_MPI(myrank,'ispec greater than nspec in mesh creation')
                endif

                ! sets index to deeper, basement layers
                ! note: doubling layer can not be at surface,
                !       since modulo 2 is zero for doubling layer, it will always have one layer below and one layer above
                !       for example: 1 doubling layer -> mesh regions = 3:
                !                    region 1 is bottom, region 2 is doubling, region 3 top
                doubling_index = IFLAG_BASEMENT_TOPO
                ispec_material_id(ispec) = material_num(ir,ix,iy)

                ! store coordinates
                call store_coords(xstore,ystore,zstore,xelm,yelm,zelm,ispec,nspec)

                ! detect mesh boundaries
                call get_flags_boundaries(nspec,iproc_xi_current,iproc_eta_current, &
                                          ispec,doubling_index, &
                                          xstore(:,:,:,ispec),ystore(:,:,:,ispec),zstore(:,:,:,ispec), &
                                          iboun,iMPIcut_xi,iMPIcut_eta,NPROC_XI,NPROC_ETA, &
                                          UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK,NEX_XI,NEX_ETA)

              enddo
            enddo
          enddo
        enddo
        ! end of loop on all the mesh points in current subregion

      endif ! regular/doubling region

    ! end of loop on all the subregions of the current region the mesh
    enddo

    ! mesh dimensions
    x_min = minval(xstore(:,:,:,:))
    x_max = maxval(xstore(:,:,:,:))
    y_min = minval(ystore(:,:,:,:))
    y_max = maxval(ystore(:,:,:,:))
    z_min = minval(zstore(:,:,:,:))
    z_max = maxval(zstore(:,:,:,:))
    call min_all_dp(x_min,x_min_glob)
    call max_all_dp(x_max,x_max_glob)
    call min_all_dp(y_min,y_min_glob)
    call max_all_dp(y_max,y_max_glob)
    call min_all_dp(z_min,z_min_glob)
    call max_all_dp(z_max,z_max_glob)
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'mesh dimensions:'
      write(IMAIN,*) '  Xmin and Xmax of the model = ',sngl(x_min_glob),sngl(x_max_glob)
      write(IMAIN,*) '  Ymin and Ymax of the model = ',sngl(y_min_glob),sngl(y_max_glob)
      write(IMAIN,*) '  Zmin and Zmax of the model = ',sngl(z_min_glob),sngl(z_max_glob)
    endif

    ! compare to exact theoretical value (bottom is always flat)
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'exact area = ',sngl((UTM_Y_MAX-UTM_Y_MIN)*(UTM_X_MAX-UTM_X_MIN)),'(m^2)'
      write(IMAIN,*) '           = ',sngl((UTM_Y_MAX-UTM_Y_MIN)*(UTM_X_MAX-UTM_X_MIN)/1000./1000.),'(km^2)'
      call flush_IMAIN()
    endif

    ! check total number of spectral elements created
    if (ispec /= nspec) call exit_MPI(myrank,'ispec should equal nspec')

    ! check that we assign a material to each element
    if (any(ispec_material_id(:) == -1000)) stop 'Element of undefined material found'

    end subroutine cmm_create_mesh_elements

!
!----------------------------------------------------------------------------------------------------
!

    subroutine cmm_create_addressing()

    implicit none
    ! local parameters
    integer :: i,j,k,ispec,ier
    integer :: ieoff,ilocnum

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'creating indirect addressing for unstructured mesh'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! puts x,y,z locations into 1D arrays
    do ispec = 1,nspec
      ieoff = NGLLCUBE_M*(ispec-1)
      ilocnum = 0
      do k = 1,NGLLZ_M
        do j = 1,NGLLY_M
          do i = 1,NGLLX_M
            ilocnum = ilocnum + 1
            xp(ilocnum+ieoff) = xstore(i,j,k,ispec)
            yp(ilocnum+ieoff) = ystore(i,j,k,ispec)
            zp(ilocnum+ieoff) = zstore(i,j,k,ispec)
          enddo
        enddo
      enddo
    enddo

    ! sorts xp,yp,zp in lexicographical order (increasing values)
    call get_global(npointot,xp,yp,zp,iglob,locval,ifseg,nglob,UTM_X_MIN,UTM_X_MAX)

    ! checks nglob range with pre-computed values
    ! note: if mesh squeezes elements such that we can't distinguish two close-by mesh points anymore
    !          the total number of mesh points might have changed
    if (nglob /= NGLOB_AB) then
      ! user output
      print *,'Error nglob: sorted value ',nglob,'differs from pre-computed ',NGLOB_AB
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Error nglob: sorted value ',nglob,'differs from pre-computed ',NGLOB_AB
        write(IMAIN,*)
        write(IMAIN,*) 'writing out problematic mesh: mesh_debug.vtk'
        write(IMAIN,*) 'please check your mesh setup...'
        call flush_IMAIN()
      endif

      ! debug file output
      ! vtk file format output
      open(IOVTK,file=prname(1:len_trim(prname))//'mesh_debug.vtk',status='unknown')
      write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
      write(IOVTK,'(a)') 'material model VTK file'
      write(IOVTK,'(a)') 'ASCII'
      write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
      write(IOVTK, '(a,i12,a)') 'POINTS ', nspec*NGLLX_M*NGLLY_M*NGLLZ_M, ' float'
      ilocnum = 0
      ibool(:,:,:,:) = 0
      do ispec = 1,nspec
        do k = 1,NGLLZ_M
          do j = 1,NGLLY_M
            do i = 1,NGLLX_M
              ilocnum = ilocnum + 1
              ibool(i,j,k,ispec) = ilocnum
              write(IOVTK,*) xstore(i,j,k,ispec),ystore(i,j,k,ispec),zstore(i,j,k,ispec)
            enddo
          enddo
        enddo
      enddo
      write(IOVTK,*)
      ! note: indices for vtk start at 0
      write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
      do ispec = 1,nspec
        write(IOVTK,'(9i12)') 8, &
              ibool(1,1,1,ispec)-1,ibool(2,1,1,ispec)-1,ibool(2,2,1,ispec)-1,ibool(1,2,1,ispec)-1, &
              ibool(1,1,2,ispec)-1,ibool(2,1,2,ispec)-1,ibool(2,2,2,ispec)-1,ibool(1,2,2,ispec)-1
      enddo
      ibool(:,:,:,:) = 0
      write(IOVTK,*)
      ! type: hexahedrons
      write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
      write(IOVTK,'(6i12)') (12,ispec=1,nspec)
      write(IOVTK,*)
      write(IOVTK,'(a,i12)') "CELL_DATA ",nspec
      write(IOVTK,'(a)') "SCALARS elem_val float"
      write(IOVTK,'(a)') "LOOKUP_TABLE default"
      do ispec = 1,nspec
        write(IOVTK,*) ispec_material_id(ispec)
      enddo
      write(IOVTK,*) ''
      close(IOVTK)

      ! stop mesher
      call exit_MPI(myrank,'incorrect global number, please check mesh input parameters')
    endif

    ! put in classical format
    allocate(nodes_coords(nglob,3),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1280')
    if (ier /= 0) stop 'Error allocating array nodes_coords'

    nodes_coords(:,:) = 0.0d0
    ibool(:,:,:,:) = 0

    do ispec = 1,nspec
      ieoff = NGLLCUBE_M*(ispec-1)
      ilocnum = 0
      do k = 1,NGLLZ_M
        do j = 1,NGLLY_M
          do i = 1,NGLLX_M
            ilocnum = ilocnum + 1
            ibool(i,j,k,ispec) = iglob(ilocnum+ieoff)
            nodes_coords(iglob(ilocnum+ieoff),1) = xstore(i,j,k,ispec)
            nodes_coords(iglob(ilocnum+ieoff),2) = ystore(i,j,k,ispec)
            nodes_coords(iglob(ilocnum+ieoff),3) = zstore(i,j,k,ispec)
          enddo
        enddo
      enddo
    enddo

    ! checks ibool range
    if (minval(ibool(:,:,:,:)) /= 1 .or. maxval(ibool(:,:,:,:)) /= nglob) then
      print *,'Error ibool: maximum value ',maxval(ibool(:,:,:,:)) ,'should be ',nglob
      call exit_MPI(myrank,'incorrect global ibool numbering')
    endif

    end subroutine cmm_create_addressing


  end subroutine create_meshfem_mesh
