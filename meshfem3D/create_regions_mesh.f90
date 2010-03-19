!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

  integer NDOUBLINGS
  integer, dimension(2) :: ner_doublings

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK
  double precision horiz_size,vert_size

  character(len=150) LOCAL_PATH

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
  character(len=150) prname,prname2

! number of elements on the boundaries
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  integer i,j,k,ia,ispec,ispec_superbrick,itype_element,ipoin
  integer iproc_xi,iproc_eta

! flag indicating whether point is in the sediments
!  logical point_is_in_sediments
  logical, dimension(:,:,:,:), allocatable :: flag_sediments
  logical, dimension(:), allocatable :: not_fully_in_bedrock

! material properties
  integer :: NMATERIALS
! first dimension  : material_id
! second dimension : #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id
  double precision , dimension(NMATERIALS,6) ::  material_properties

! doublings zone
  integer :: nspec_sb
  integer :: offset_x,offset_y,offset_z
  logical, dimension(NSPEC_DOUBLING_SUPERBRICK,6) :: iboun_sb
  integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick
  double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick


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



!--- apply heuristic rule to modify doubling regions to balance angles


!   if(APPLY_HEURISTIC_RULE .and. .not. USE_REGULAR_MESH) then

!      stop 'pas encore implemente'

! ! define number of subregions affected by heuristic rule in doubling regions
!   nsubregions = 8

!   do isubregion = 1,nsubregions

! ! define shape of elements for heuristic
!     call define_subregions_heuristic(myrank,isubregion,iaddx,iaddy,iaddz, &
!          nblayers,ner_layer, &
!          ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
!          itype_element,npx,npy)
!          !              NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM)

! ! loop on all the mesh points in current subregion
!   do ir = ir1,ir2,dir
!     do iy = iy1,iy2,diy
!       do ix = ix1,ix2,dix

! ! this heuristic rule is only valid for 8-node elements
! ! it would not work in the case of 27 nodes

! !----
!     if(itype_element == ITYPE_UNUSUAL_1) then

! ! side 1
!       horiz_size = xgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) &
!                  - xgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
!       xgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) = &
!          xgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + horiz_size * MAGIC_RATIO

!       vert_size = zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) &
!                  - zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
!       zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) = &
!          zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + vert_size * MAGIC_RATIO / 0.50

! ! side 2
!       horiz_size = xgrid(ir+iar*iaddz(3),ix+iax*iaddx(3),iy+iay*iaddy(3)) &
!                  - xgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4))
!       xgrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) = &
!          xgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) + horiz_size * MAGIC_RATIO

!       vert_size = zgrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) &
!                  - zgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4))
!       zgrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) = &
!          zgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) + vert_size * MAGIC_RATIO / 0.50

! !----
!     else if(itype_element == ITYPE_UNUSUAL_1p) then

! ! side 1
!       horiz_size = xgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) &
!                  - xgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
!       xgrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) = &
!          xgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + horiz_size * (1. - MAGIC_RATIO)

!       vert_size = zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) &
!                  - zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
!       zgrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) = &
!          zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + vert_size * MAGIC_RATIO / 0.50

! ! side 2
!       horiz_size = xgrid(ir+iar*iaddz(3),ix+iax*iaddx(3),iy+iay*iaddy(3)) &
!                  - xgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4))
!       xgrid(ir+iar*iaddz(7),ix+iax*iaddx(7),iy+iay*iaddy(7)) = &
!          xgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) + horiz_size * (1. - MAGIC_RATIO)

!       vert_size = zgrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) &
!                  - zgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4))
!       zgrid(ir+iar*iaddz(7),ix+iax*iaddx(7),iy+iay*iaddy(7)) = &
!          zgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) + vert_size * MAGIC_RATIO / 0.50

! !----
!     else if(itype_element == ITYPE_UNUSUAL_4) then

! ! side 1
!       horiz_size = ygrid(ir+iar*iaddz(3),ix+iax*iaddx(3),iy+iay*iaddy(3)) &
!                  - ygrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2))
!       ygrid(ir+iar*iaddz(7),ix+iax*iaddx(7),iy+iay*iaddy(7)) = &
!          ygrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) + horiz_size * (1. - MAGIC_RATIO)

!       vert_size = zgrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) &
!                  - zgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2))
!       zgrid(ir+iar*iaddz(7),ix+iax*iaddx(7),iy+iay*iaddy(7)) = &
!          zgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) + vert_size * MAGIC_RATIO / 0.50

! ! side 2
!       horiz_size = ygrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) &
!                  - ygrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
!       ygrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) = &
!          ygrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + horiz_size * (1. - MAGIC_RATIO)

!       vert_size = zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) &
!                  - zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
!       zgrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) = &
!          zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + vert_size * MAGIC_RATIO / 0.50

! !----
!     else if(itype_element == ITYPE_UNUSUAL_4p) then

! ! side 1
!       horiz_size = ygrid(ir+iar*iaddz(3),ix+iax*iaddx(3),iy+iay*iaddy(3)) &
!                  - ygrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2))
!       ygrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) = &
!          ygrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) + horiz_size * MAGIC_RATIO

!       vert_size = zgrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) &
!                  - zgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2))
!       zgrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) = &
!          zgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) + vert_size * MAGIC_RATIO / 0.50

! ! side 2
!       horiz_size = ygrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) &
!                  - ygrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
!       ygrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) = &
!          ygrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + horiz_size * MAGIC_RATIO

!       vert_size = zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) &
!                  - zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
!       zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) = &
!          zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + vert_size * MAGIC_RATIO / 0.50

!    endif

! enddo
! enddo
! enddo

! enddo

! endif

!---

! generate the elements in all the regions of the mesh
  ispec = 0

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

! A revoir
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

               ! A revoir
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
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  write(prname2, "(i6.6,'.dx')") myrank
  open(unit=66,file='DX_HR_'//prname2,status='unknown')
  
! ************* generate points ******************

!   NPOIN = 0
! ! loop on all the mesh points
!   do ir = 0,2*NER!ir1,ir2,dir
!      do iy = 0,2*NEX_PER_PROC_ETA!iy1,iy2,diy
!         do ix = 0,2*NEX_PER_PROC_XI!ix1,ix2,dix
!            !if(grid_flag_MR(ir,ix,iy)) then
!               NPOIN = NPOIN + 1
!            !end if
!         end do
!      end do
!   end do

! ! write OpenDX header
!   write(66,*) 'object 1 class array type float rank 1 shape 3 items ',NPOIN,' data follows'

! !! 0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA
! ! loop on all the mesh points
!   do ir = 0,2*NER!ir1,ir2,dir
!      do iy = 0,2*NEX_PER_PROC_ETA!iy1,iy2,diy
!         do ix = 0,2*NEX_PER_PROC_XI!ix1,ix2,dix
!            !if(grid_flag_MR(ir,ix,iy)) then
!               write(66,*) sngl(xgrid(ir,ix,iy)),sngl(ygrid(ir,ix,iy)),sngl(zgrid(ir,ix,iy))
!            !end if
!         end do
!      end do
!   end do


! write OpenDX header
  write(66,*) 'object 1 class array type float rank 1 shape 3 items ',nglob,' data follows'
  
  do ipoin = 1,nglob
     write(66,*) sngl(nodes_coords(ipoin,1)),sngl(nodes_coords(ipoin,2)),sngl(nodes_coords(ipoin,3))
  end do

! ************* generate elements ******************

  write(66,*) 'object 2 class array type int rank 1 shape ',8,' items ',nspec,' data follows'

  do ispec=1,nspec

        ! point order in OpenDX in 2D is 1,4,2,3 *not* 1,2,3,4 as in AVS
        ! point order in OpenDX in 3D is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
        ! in the case of OpenDX, node numbers start at zero     
        write(66,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") &
             ibool(1,1,2,ispec)-1,ibool(2,1,2,ispec)-1,ibool(1,2,2,ispec)-1,ibool(2,2,2,ispec)-1,&
             ibool(1,1,1,ispec)-1,ibool(2,1,1,ispec)-1,ibool(1,2,1,ispec)-1,ibool(2,2,1,ispec)-1
!              ibool(5,1,5,ispec)-1,ibool(5,1,1,ispec)-1,ibool(1,1,1,ispec)-1,ibool(1,1,5,ispec)-1,&
!              ibool(5,5,5,ispec)-1,ibool(5,5,1,ispec)-1,ibool(1,5,1,ispec)-1,ibool(1,5,5,ispec)-1
!             ibool(1,5,1,ispec)-1,ibool(1,1,1,ispec)-1,ibool(1,5,5,ispec)-1,ibool(1,1,5,ispec)-1,&
!             ibool(5,5,1,ispec)-1,ibool(5,1,1,ispec)-1,ibool(5,5,5,ispec)-1,ibool(5,1,5,ispec)-1
!             ibool(1,5,5,ispec)-1,ibool(1,1,5,ispec)-1,ibool(5,5,5,ispec)-1,ibool(5,1,5,ispec)-1,&
!             ibool(1,5,1,ispec)-1,ibool(1,1,1,ispec)-1,ibool(5,5,1,ispec)-1,ibool(5,1,1,ispec)-1
!             ibool(1,1,1,ispec)-1,ibool(1,5,1,ispec)-1,ibool(5,5,1,ispec)-1,ibool(5,1,1,ispec)-1,&
!             ibool(1,1,5,ispec)-1,ibool(1,5,5,ispec)-1,ibool(5,5,5,ispec)-1,ibool(5,1,5,ispec)-1
        !             point 1 = (0,0,0), point 2 = (0,1,0), point 3 = (1,1,0), point 4 = (1,0,0)
        !          then top (positive z-direction) of element
        !             point 5 = (0,0,1), point 6 = (0,1,1), point 7 = (1,1,1), point 8 = (1,0,1)
        !          ibool(4,ispec)-1, ibool(1,ispec)-1, ibool(8,ispec)-1, ibool(5,ispec)-1, &
        !          ibool(3,ispec)-1, ibool(2,ispec)-1, ibool(7,ispec)-1, ibool(6,ispec)-1      
     
  enddo



! ************* generate element data values ******************

! output OpenDX header for data
  write(66,*) 'attribute "element type" string "cubes"'

  write(66,*) 'attribute "ref" string "positions"'
  write(66,*) 'object 3 class array type float rank 0 items ',nspec,' data follows'

! loop on all the elements
  do ispec = 1,nspec
    write(66,*) 1.
  enddo


  write(66,*) 'attribute "dep" string "connections"'
  write(66,*) 'object "irregular positions irregular connections" class field'
  write(66,*) 'component "positions" value 1'
  write(66,*) 'component "connections" value 2'
  write(66,*) 'component "data" value 3'
  write(66,*) 'end'

  close(66)
  call sync_all()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(minval(ibool(:,:,:,:)) /= 1 .or. maxval(ibool(:,:,:,:)) /= NGLOB_AB) then
     print*,maxval(ibool(:,:,:,:)) ,NGLOB_AB
    call exit_MPI(myrank,'incorrect global numbering')
 end if



  call store_boundaries(myrank,iboun,nspec, &
      ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
      nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
      NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

  call save_databases(prname,nspec,nglob,iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta,&
       ibool,nodes_coords,true_material_num,nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,&
       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,&
       NMATERIALS,material_properties)

  end subroutine create_regions_mesh

end module createRegMesh


