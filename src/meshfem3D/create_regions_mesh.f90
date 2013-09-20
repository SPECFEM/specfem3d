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

module createRegMesh
contains

  subroutine create_regions_mesh(xgrid,ygrid,zgrid,ibool, &
                               xstore,ystore,zstore,iproc_xi,iproc_eta,addressing,nspec, &
                               NGLOB_AB,npointot, &
                               NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER, &
                               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                               NPROC_XI,NPROC_ETA, &
                               nsubregions,subregions,nblayers,ner_layer,NMATERIALS,material_properties, &
                               myrank, sizeprocs, &
                 LOCAL_PATH,UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
                               CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES, &
                               USE_REGULAR_MESH,NDOUBLINGS,ner_doublings, &
                               ADIOS_ENABLED, ADIOS_FOR_DATABASES)

    ! create the different regions of the mesh
  use adios_manager_mod
  use mpi

    implicit none

    include "constants.h"
    include "constants_meshfem3D.h"

    ! number of spectral elements in each block
    integer nspec

    integer NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER

    integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP

    integer NPROC_XI,NPROC_ETA

    integer npointot

    logical USE_REGULAR_MESH
    logical CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES
    logical ADIOS_ENABLED, ADIOS_FOR_DATABASES

    integer NDOUBLINGS
    integer, dimension(2) :: ner_doublings

    double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK
    !double precision horiz_size,vert_size

    integer addressing(0:NPROC_XI-1,0:NPROC_ETA-1)

    ! arrays with the mesh
    double precision xgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA)
    double precision ygrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA)
    double precision zgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA)

    double precision, dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec) :: xstore,ystore,zstore

    integer,dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec) :: ibool

    character(len=256) :: LOCAL_PATH

    ! auxiliary variables to generate the mesh
    !  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: nodes_coords
    double precision, dimension(:,:), allocatable :: nodes_coords
    integer ix,iy,ir,ir1,ir2,dir
    integer ix1,ix2,dix,iy1,iy2,diy
    integer iax,iay,iar
    integer isubregion,nsubregions,doubling_index,nmeshregions
    integer imaterial_number
    integer, dimension(nspec) :: true_material_num
    integer, dimension(:,:,:), allocatable :: material_num
    !integer material_num(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA)

    !  definition of the different regions of the model in the mesh (nx,ny,nz)
    !  #1 #2 : nx_begining,nx_end
    !  #3 #4 : ny_begining,ny_end
    !  #5 #6 : nz_begining,nz_end
    !     #7 : material number
    integer subregions(nsubregions,7)

    ! layers of the model
    integer nblayers
    integer ner_layer(nblayers)

    ! topology of the elements
    integer iaddx(NGNOD_EIGHT_CORNERS)
    integer iaddy(NGNOD_EIGHT_CORNERS)
    integer iaddz(NGNOD_EIGHT_CORNERS)

    double precision xelm(NGNOD_EIGHT_CORNERS)
    double precision yelm(NGNOD_EIGHT_CORNERS)
    double precision zelm(NGNOD_EIGHT_CORNERS)

    ! boundary locator
    logical, dimension(:,:), allocatable :: iboun

    ! proc numbers for MPI
    integer myrank, sizeprocs

    ! variables for creating array ibool (some arrays also used for AVS or DX files)
    integer, dimension(:), allocatable :: iglob,locval
    logical, dimension(:), allocatable :: ifseg
    double precision, dimension(:), allocatable :: xp,yp,zp

    integer nglob,NGLOB_AB
    integer ieoff,ilocnum

    ! boundary parameters locator
    integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top

    ! MPI cut-planes parameters along xi and along eta
    logical, dimension(:,:), allocatable :: iMPIcut_xi,iMPIcut_eta

    ! name of the database file
    character(len=256) prname !,prname2

    ! number of elements on the boundaries
    integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

    integer i,j,k,ia,ispec,ispec_superbrick,ier ! itype_element ,ipoin
    integer iproc_xi,iproc_eta

    ! flag indicating whether point is in the sediments
    !  logical point_is_in_sediments
    logical, dimension(:,:,:,:), allocatable :: flag_sediments
    logical, dimension(:), allocatable :: not_fully_in_bedrock

    ! material properties
    integer :: NMATERIALS
    double precision :: VP_MAX
    ! first dimension  : material_id
    ! second dimension : #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id
    double precision , dimension(NMATERIALS,6) ::  material_properties

    ! doublings zone
    integer :: nspec_sb
    integer :: ioffset_x,ioffset_y,ioffset_z
    logical, dimension(NSPEC_DOUBLING_SUPERBRICK,6) :: iboun_sb
    integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick
    double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick

    logical, parameter :: DEBUG = .true.

    ! **************

    ! create the name for the database of the current slide and region
    call create_name_database(prname,myrank,LOCAL_PATH)

    ! flag indicating whether point is in the sediments
    allocate(flag_sediments(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array flag_sediments'
    allocate(not_fully_in_bedrock(nspec),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'not enough memory to allocate arrays')

    ! boundary locator
    allocate(iboun(6,nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array iboun'

    ! boundary parameters locator
    allocate(ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_xmin'
    allocate(ibelm_xmax(NSPEC2DMAX_XMIN_XMAX),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_xmax'
    allocate(ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_ymin'
    allocate(ibelm_ymax(NSPEC2DMAX_YMIN_YMAX),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_ymax'
    allocate(ibelm_bottom(NSPEC2D_BOTTOM),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibelm_bottom'
    allocate(ibelm_top(NSPEC2D_TOP),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'not enough memory to allocate arrays')

    ! MPI cut-planes parameters along xi and along eta
    allocate(iMPIcut_xi(2,nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating array iMPIcut_xi'
    allocate(iMPIcut_eta(2,nspec),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'not enough memory to allocate arrays')

    ! allocate memory for arrays
    allocate(iglob(npointot),stat=ier)
    if( ier /= 0 ) stop 'error allocating array iglob'
    allocate(locval(npointot),stat=ier)
    if( ier /= 0 ) stop 'error allocating array locval'
    allocate(ifseg(npointot),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ifseg'
    allocate(xp(npointot),stat=ier)
    if( ier /= 0 ) stop 'error allocating array xp'
    allocate(yp(npointot),stat=ier)
    if( ier /= 0 ) stop 'error allocating array yp'
    allocate(zp(npointot),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'not enough memory to allocate arrays')

    ! allocate material ids array
    allocate(material_num(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'not enough memory to allocate arrays')

    ! generate the elements in all the regions of the mesh
    ispec = 0

    ! to check that we assign a material to each element
    material_num(:,:,:) = -1000 ! dummy value

    do isubregion = 1,nsubregions
       call define_model_regions(NEX_PER_PROC_XI,NEX_PER_PROC_ETA,iproc_xi,iproc_eta,&
            isubregion,nsubregions,subregions,nblayers,ner_layer, &
            iaddx,iaddy,iaddz,ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
            imaterial_number)

       ! loop on all the mesh points in current subregion
       do ir = ir1,ir2,dir
          do iy = iy1,iy2,diy
             do ix = ix1,ix2,dix

                material_num(ir,ix,iy) = imaterial_number

                ! end of loop on all the mesh points in current subregion
             enddo
          enddo
       enddo

       ! end of loop on all the subregions of the current region the mesh
    enddo

    if(USE_REGULAR_MESH) then
       nmeshregions = 1
    else
       call define_superbrick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb)
       nspec_sb = NSPEC_DOUBLING_SUPERBRICK
       if(NDOUBLINGS == 1) then
          nmeshregions = 3
       else if(NDOUBLINGS == 2) then
          nmeshregions = 5
       else
          stop 'Wrong number of mesh regions'
       endif
    endif

    do isubregion = 1,nmeshregions
      ! define shape of elements
      call define_mesh_regions(USE_REGULAR_MESH,isubregion,NER,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,iproc_xi,iproc_eta,&
            nblayers,ner_layer,NDOUBLINGS,ner_doublings,&
            iaddx,iaddy,iaddz,ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar)

      ! loop on all the mesh points in current subregion
      do ir = ir1,ir2,dir
        do iy = iy1,iy2,diy
          do ix = ix1,ix2,dix

            if(modulo(isubregion,2) == 1) then

              ! Regular subregion case

              ! loop over the NGNOD_EIGHT_CORNERS nodes
              do ia=1,NGNOD_EIGHT_CORNERS
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
              if(ispec > nspec) then
                call exit_MPI(myrank,'ispec greater than nspec in mesh creation')
              endif

              ! check if element is on topography
              if((ir == ir2) .and. (isubregion == nmeshregions)) then
                doubling_index = IFLAG_ONE_LAYER_TOPOGRAPHY
              else
                doubling_index = IFLAG_BASEMENT_TOPO
              endif

              true_material_num(ispec) = material_num(ir,ix,iy)

              ! store coordinates
              call store_coords(xstore,ystore,zstore,xelm,yelm,zelm,ispec,nspec)

              ! detect mesh boundaries
              call get_flags_boundaries(nspec,iproc_xi,iproc_eta,ispec,doubling_index, &
                  xstore(:,:,:,ispec),ystore(:,:,:,ispec),zstore(:,:,:,ispec), &
                  iboun,iMPIcut_xi,iMPIcut_eta,NPROC_XI,NPROC_ETA, &
                  UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK)

            else

              ! Irregular subregion case

              ! loop on all the elements in the mesh doubling superbrick
              do ispec_superbrick = 1,nspec_sb
                ! loop on all the corner nodes of this element
                do ia = 1,NGNOD_EIGHT_CORNERS
                  ! define topological coordinates of this mesh point
                  ioffset_x = ix + iax*x_superbrick(ibool_superbrick(ia,ispec_superbrick))
                  ioffset_y = iy + iay*y_superbrick(ibool_superbrick(ia,ispec_superbrick))
                  ioffset_z = ir + iar*z_superbrick(ibool_superbrick(ia,ispec_superbrick))

                  xelm(ia) = xgrid(ioffset_z,ioffset_x,ioffset_y)
                  yelm(ia) = ygrid(ioffset_z,ioffset_x,ioffset_y)
                  zelm(ia) = zgrid(ioffset_z,ioffset_x,ioffset_y)

                enddo

                ! add one spectral element to the list and store its material number
                ispec = ispec + 1
                if(ispec > nspec) then
                  call exit_MPI(myrank,'ispec greater than nspec in mesh creation')
                endif

                doubling_index = IFLAG_BASEMENT_TOPO
                true_material_num(ispec) = material_num(ir,ix,iy)

                ! store coordinates
                call store_coords(xstore,ystore,zstore,xelm,yelm,zelm,ispec,nspec)

                ! detect mesh boundaries
                call get_flags_boundaries(nspec,iproc_xi,iproc_eta,ispec,doubling_index, &
                     xstore(:,:,:,ispec),ystore(:,:,:,ispec),zstore(:,:,:,ispec), &
                     iboun,iMPIcut_xi,iMPIcut_eta,NPROC_XI,NPROC_ETA, &
                     UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK)

              enddo
            endif

          ! end of loop on all the mesh points in current subregion
          enddo
        enddo
      enddo

    ! end of loop on all the subregions of the current region the mesh
    enddo

    ! check total number of spectral elements created
    if(ispec /= nspec) call exit_MPI(myrank,'ispec should equal nspec')

    ! check that we assign a material to each element
    if(any(true_material_num(:) == -1000)) stop 'Element of undefined material found'

    ! puts x,y,z locations into 1D arrays
    do ispec=1,nspec
       ieoff = NGLLCUBE_M*(ispec-1)
       ilocnum = 0
       do k=1,NGLLZ_M
          do j=1,NGLLY_M
             do i=1,NGLLX_M
                ilocnum = ilocnum + 1
                xp(ilocnum+ieoff) = xstore(i,j,k,ispec)
                yp(ilocnum+ieoff) = ystore(i,j,k,ispec)
                zp(ilocnum+ieoff) = zstore(i,j,k,ispec)
             enddo
          enddo
       enddo
    enddo

    ! sorts xp,yp,zp in lexicographical order (increasing values)
    call get_global(nspec,xp,yp,zp,iglob,locval,ifseg,nglob,npointot,UTM_X_MIN,UTM_X_MAX)

    ! checks nglob range with pre-computed values
    ! note: if mesh squeezes elements such that we can't distinguish two close-by mesh points anymore
    !          the total number of mesh points might have changed
    if(nglob /= NGLOB_AB) then
      ! user output
      print*,'error nglob: sorted value ',nglob,'differs from pre-computed ',NGLOB_AB
      if(myrank == 0 ) then
        write(IMAIN,*)
        write(IMAIN,*) 'error nglob: sorted value ',nglob,'differs from pre-computed ',NGLOB_AB
        write(IMAIN,*)
        write(IMAIN,*) 'writing out problematic mesh: mesh_debug.vtk'
        write(IMAIN,*) 'please check your mesh setup...'
      endif
      ! debug file output
      ! vtk file format output
      open(66,file=prname(1:len_trim(prname))//'mesh_debug.vtk',status='unknown')
      write(66,'(a)') '# vtk DataFile Version 3.1'
      write(66,'(a)') 'material model VTK file'
      write(66,'(a)') 'ASCII'
      write(66,'(a)') 'DATASET UNSTRUCTURED_GRID'
      write(66, '(a,i12,a)') 'POINTS ', nspec*NGLLX_M*NGLLY_M*NGLLZ_M, ' float'
      ilocnum = 0
      ibool(:,:,:,:) = 0
      do ispec=1,nspec
        do k=1,NGLLZ_M
          do j=1,NGLLY_M
            do i=1,NGLLX_M
              ilocnum = ilocnum + 1
              ibool(i,j,k,ispec) = ilocnum
              write(66,*) xstore(i,j,k,ispec),ystore(i,j,k,ispec),zstore(i,j,k,ispec)
            enddo
          enddo
        enddo
      enddo
      write(66,*)
      ! note: indices for vtk start at 0
      write(66,'(a,i12,i12)') "CELLS ",nspec,nspec*9
      do ispec=1,nspec
        write(66,'(9i12)') 8, &
              ibool(1,1,1,ispec)-1,ibool(2,1,1,ispec)-1,ibool(2,2,1,ispec)-1,ibool(1,2,1,ispec)-1,&
              ibool(1,1,2,ispec)-1,ibool(2,1,2,ispec)-1,ibool(2,2,2,ispec)-1,ibool(1,2,2,ispec)-1
      enddo
      ibool(:,:,:,:) = 0
      write(66,*)
      ! type: hexahedrons
      write(66,'(a,i12)') "CELL_TYPES ",nspec
      write(66,'(6i12)') (12,ispec=1,nspec)
      write(66,*)
      write(66,'(a,i12)') "CELL_DATA ",nspec
      write(66,'(a)') "SCALARS elem_val float"
      write(66,'(a)') "LOOKUP_TABLE default"
      do ispec = 1,nspec
        write(66,*) true_material_num(ispec)
      enddo
      write(66,*) ""
      close(66)

      ! stop mesher
      call exit_MPI(myrank,'incorrect global number, please check mesh input parameters')
    endif

    ! put in classical format
    allocate(nodes_coords(nglob,3),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_coords'
    nodes_coords(:,:) = 0.0d0
    ibool(:,:,:,:) = 0

    do ispec=1,nspec
       ieoff = NGLLCUBE_M*(ispec-1)
       ilocnum = 0
       do k=1,NGLLZ_M
          do j=1,NGLLY_M
             do i=1,NGLLX_M
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
    if(minval(ibool(:,:,:,:)) /= 1 .or. maxval(ibool(:,:,:,:)) /= nglob) then
       print*,'error ibool: maximum value ',maxval(ibool(:,:,:,:)) ,'should be ',nglob
       call exit_MPI(myrank,'incorrect global ibool numbering')
    endif

    !--- Initialize ADIOS and setup the buffer size
    if (ADIOS_ENABLED) then
      call adios_setup()
    endif

    ! outputs mesh file for visualization
    call create_visual_files(CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES, &
                            nspec,nglob, &
                            prname,nodes_coords,ibool,true_material_num)

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
       nspec,nglob,iproc_xi,iproc_eta, &
       NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta,&
       ibool,nodes_coords,true_material_num, &
       nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
       NSPEC2D_BOTTOM,NSPEC2D_TOP, NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
       ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top,&
       NMATERIALS,material_properties)
  else
    ! saves mesh as databases file
    call save_databases(prname,nspec,nglob,iproc_xi,iproc_eta, &
                      NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta,&
                      ibool,nodes_coords,true_material_num, &
                      nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP,&
                      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                      ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top,&
                      NMATERIALS,material_properties)
  endif

    !--- Clean ADIOS. Make sure everything is already written
    if (ADIOS_ENABLED) then
      call adios_cleanup()
    endif

  end subroutine create_regions_mesh

end module createRegMesh

