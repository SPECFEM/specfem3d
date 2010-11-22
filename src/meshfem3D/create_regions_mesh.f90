!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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



module createRegMesh
contains


  subroutine create_regions_mesh(xgrid,ygrid,zgrid,ibool, &
       xstore,ystore,zstore,npx,npy,iproc_xi,iproc_eta,addressing,nspec, &
       NGLOB_AB,npointot, &
       NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER, &
       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
       NPROC_XI,NPROC_ETA,NSPEC2D_A_XI,NSPEC2D_B_XI, &
       NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
       nsubregions,subregions,nblayers,ner_layer,NMATERIALS,material_properties, &
       myrank,LOCAL_PATH,UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
       CREATE_ABAQUS_FILES,CREATE_DX_FILES,&
       USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

    ! create the different regions of the mesh

    implicit none

    include "constants.h"

    ! number of spectral elements in each block
    integer nspec

    integer NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER

    integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP

    integer NPROC_XI,NPROC_ETA,NSPEC2D_A_XI,NSPEC2D_B_XI
    integer NSPEC2D_A_ETA,NSPEC2D_B_ETA

    integer npx,npy
    integer npointot

    logical USE_REGULAR_MESH
    logical CREATE_ABAQUS_FILES,CREATE_DX_FILES

    integer NDOUBLINGS
    integer, dimension(2) :: ner_doublings

    double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK
    !double precision horiz_size,vert_size

    character(len=256) LOCAL_PATH

    integer addressing(0:NPROC_XI-1,0:NPROC_ETA-1)

    ! arrays with the mesh
    double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore
    !  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: nodes_coords
    double precision, dimension(:,:), allocatable :: nodes_coords

    double precision xgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA)
    double precision ygrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA)
    double precision zgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA)

    integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

    ! auxiliary variables to generate the mesh
    integer ix,iy,ir,ir1,ir2,dir
    integer ix1,ix2,dix,iy1,iy2,diy
    integer iax,iay,iar
    integer isubregion,nsubregions,doubling_index,nmeshregions
    integer imaterial_number
    integer true_material_num(nspec)
    integer material_num(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA)

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
    integer iaddx(NGNOD)
    integer iaddy(NGNOD)
    integer iaddz(NGNOD)

    double precision xelm(NGNOD)
    double precision yelm(NGNOD)
    double precision zelm(NGNOD)

    ! boundary locator
    logical, dimension(:,:), allocatable :: iboun

    ! proc numbers for MPI
    integer myrank

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

    integer i,j,k,ia,ispec,ispec_superbrick ! itype_element ,ipoin
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
    integer :: offset_x,offset_y,offset_z
    logical, dimension(NSPEC_DOUBLING_SUPERBRICK,6) :: iboun_sb
    integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick
    double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick

    ! to avoid compiler warnings
    integer idummy
    idummy = NSPEC2D_A_ETA
    idummy = NSPEC2D_B_ETA
    idummy = NSPEC2D_A_XI
    idummy = NSPEC2D_B_XI
    idummy = npx
    idummy = npy

    ! **************

    ! create the name for the database of the current slide and region
    call create_name_database(prname,myrank,LOCAL_PATH)

    ! flag indicating whether point is in the sediments
    allocate(flag_sediments(NGLLX,NGLLY,NGLLZ,nspec))
    allocate(not_fully_in_bedrock(nspec))

    ! boundary locator
    allocate(iboun(6,nspec))

    ! boundary parameters locator
    allocate(ibelm_xmin(NSPEC2DMAX_XMIN_XMAX))
    allocate(ibelm_xmax(NSPEC2DMAX_XMIN_XMAX))
    allocate(ibelm_ymin(NSPEC2DMAX_YMIN_YMAX))
    allocate(ibelm_ymax(NSPEC2DMAX_YMIN_YMAX))
    allocate(ibelm_bottom(NSPEC2D_BOTTOM))
    allocate(ibelm_top(NSPEC2D_TOP))

    ! MPI cut-planes parameters along xi and along eta
    allocate(iMPIcut_xi(2,nspec))
    allocate(iMPIcut_eta(2,nspec))

    ! allocate memory for arrays
    allocate(iglob(npointot))
    allocate(locval(npointot))
    allocate(ifseg(npointot))
    allocate(xp(npointot))
    allocate(yp(npointot))
    allocate(zp(npointot))

    ! generate the elements in all the regions of the mesh
    ispec = 0


    ! to check that we assign a material to each element
    material_num(:,:,:) = -1000

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
       !call define_superbrick_one_layer(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb)
       nspec_sb = NSPEC_DOUBLING_SUPERBRICK
       !nspec_sb = NSPEC_SUPERBRICK_1L
       if(NDOUBLINGS == 1) then
          nmeshregions = 3
       else if(NDOUBLINGS == 2) then
          nmeshregions = 5
       else
          stop 'Wrong number of mesh regions'
       end if
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

                if(modulo(isubregion,2) == 1) then    ! Regular subregion case

                   ! loop over the NGNOD nodes
                   do ia=1,NGNOD
                      xelm(ia) = xgrid(ir+iar*iaddz(ia),ix+iax*iaddx(ia),iy+iay*iaddy(ia))
                      yelm(ia) = ygrid(ir+iar*iaddz(ia),ix+iax*iaddx(ia),iy+iay*iaddy(ia))
                      zelm(ia) = zgrid(ir+iar*iaddz(ia),ix+iax*iaddx(ia),iy+iay*iaddy(ia))
                   enddo

                   ! add one spectral element to the list and store its material number
                   ispec = ispec + 1
                   if(ispec > nspec) then
                      call exit_MPI(myrank,'ispec greater than nspec in mesh creation')
                   end if

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

                else             ! Irregular subregion case

                   ! loop on all the elements in the mesh doubling superbrick
                   do ispec_superbrick = 1,nspec_sb
                      ! loop on all the corner nodes of this element
                      do ia = 1,NGNOD_EIGHT_CORNERS
                         ! define topological coordinates of this mesh point
                         offset_x = ix + iax*x_superbrick(ibool_superbrick(ia,ispec_superbrick))
                         offset_y = iy + iay*y_superbrick(ibool_superbrick(ia,ispec_superbrick))
                         offset_z = ir + iar*z_superbrick(ibool_superbrick(ia,ispec_superbrick))

                         xelm(ia) = xgrid(offset_z,offset_x,offset_y)
                         yelm(ia) = ygrid(offset_z,offset_x,offset_y)
                         zelm(ia) = zgrid(offset_z,offset_x,offset_y)
                      enddo

                      ! add one spectral element to the list and store its material number
                      ispec = ispec + 1
                      if(ispec > nspec) then
                         call exit_MPI(myrank,'ispec greater than nspec in mesh creation')
                      end if


                      doubling_index = IFLAG_BASEMENT_TOPO
                      true_material_num(ispec) = material_num(ir,ix,iy)

                      ! store coordinates
                      call store_coords(xstore,ystore,zstore,xelm,yelm,zelm,ispec,nspec)

                      ! detect mesh boundaries
                      call get_flags_boundaries(nspec,iproc_xi,iproc_eta,ispec,doubling_index, &
                           xstore(:,:,:,ispec),ystore(:,:,:,ispec),zstore(:,:,:,ispec), &
                           iboun,iMPIcut_xi,iMPIcut_eta,NPROC_XI,NPROC_ETA, &
                           UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK)

                   end do
                end if

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

    do ispec=1,nspec
       ieoff = NGLLCUBE*(ispec-1)
       ilocnum = 0
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                ilocnum = ilocnum + 1
                xp(ilocnum+ieoff) = xstore(i,j,k,ispec)
                yp(ilocnum+ieoff) = ystore(i,j,k,ispec)
                zp(ilocnum+ieoff) = zstore(i,j,k,ispec)
             enddo
          enddo
       enddo
    enddo

    call get_global(nspec,xp,yp,zp,iglob,locval,ifseg,nglob,npointot,UTM_X_MIN,UTM_X_MAX)

    ! put in classical format
    allocate(nodes_coords(nglob,3))

    do ispec=1,nspec
       ieoff = NGLLCUBE*(ispec-1)
       ilocnum = 0
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                ilocnum = ilocnum + 1
                ibool(i,j,k,ispec) = iglob(ilocnum+ieoff)
                nodes_coords(iglob(ilocnum+ieoff),1) = xstore(i,j,k,ispec)
                nodes_coords(iglob(ilocnum+ieoff),2) = ystore(i,j,k,ispec)
                nodes_coords(iglob(ilocnum+ieoff),3) = zstore(i,j,k,ispec)
             enddo
          enddo
       enddo
    enddo

    if(minval(ibool(:,:,:,:)) /= 1 .or. maxval(ibool(:,:,:,:)) /= NGLOB_AB) then
       print*,maxval(ibool(:,:,:,:)) ,NGLOB_AB
       call exit_MPI(myrank,'incorrect global numbering')
    end if

    call create_visual_files(CREATE_ABAQUS_FILES,CREATE_DX_FILES,nspec,nglob,prname,nodes_coords,ibool,true_material_num)

    call store_boundaries(myrank,iboun,nspec, &
         ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
         nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
         NSPEC2D_BOTTOM,NSPEC2D_TOP, &
         NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

    VP_MAX = maxval(material_properties(:,2))
    call check_mesh_quality(myrank,VP_MAX,nglob,nspec,nodes_coords(:,1),nodes_coords(:,2),nodes_coords(:,3),ibool)

      call save_databases(prname,nspec,nglob,iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta,&
           ibool,nodes_coords,true_material_num,nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP,&
           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top,&
           NMATERIALS,material_properties)

  end subroutine create_regions_mesh

end module createRegMesh


