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

  subroutine define_model_regions(NEX_PER_PROC_XI,NEX_PER_PROC_ETA,iproc_xi,iproc_eta,&
       isubregion,nbsubregions,subregions,nblayers,ner_layer,&
       iaddx,iaddy,iaddz,ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
       num_material)

    implicit none

    include "constants.h"

    integer NEX_PER_PROC_XI,NEX_PER_PROC_ETA
    integer iproc_xi,iproc_eta

    integer ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir
    integer iax,iay,iar
    integer isubregion,nbsubregions,nblayers
    integer num_material

! topology of the elements
    integer iaddx(NGNOD_EIGHT_CORNERS)
    integer iaddy(NGNOD_EIGHT_CORNERS)
    integer iaddz(NGNOD_EIGHT_CORNERS)
    integer ner_layer(nblayers)

!  definition of the different regions of the model in the mesh (nx,ny,nz)
!  #1 #2 : nx_begining,nx_end
!  #3 #4 : ny_begining,ny_end
!  #5 #6 : nz_begining,nz_end
!     #7 : material number
    integer subregions(nbsubregions,7)

    ! to avoid compiler warnings
    integer idummy
    idummy = ner_layer(1)
! **************

     call usual_hex_nodes(NGNOD_EIGHT_CORNERS,iaddx,iaddy,iaddz)


     ix1=2*(subregions(isubregion,1) - iproc_xi*NEX_PER_PROC_XI - 1)
     if(ix1 < 0) ix1 = 0
     ix2=2*(subregions(isubregion,2) - iproc_xi*NEX_PER_PROC_XI - 1)
     if(ix2 > 2*(NEX_PER_PROC_XI - 1)) ix2 = 2*(NEX_PER_PROC_XI - 1)
     dix=2

     iy1=2*(subregions(isubregion,3) - iproc_eta*NEX_PER_PROC_ETA - 1)
     if(iy1 < 0) iy1 = 0
     iy2=2*(subregions(isubregion,4) - iproc_eta*NEX_PER_PROC_ETA - 1)
     if(iy2 > 2*(NEX_PER_PROC_ETA - 1)) iy2 = 2*(NEX_PER_PROC_ETA - 1)
     diy=2

     ir1=2*(subregions(isubregion,5) - 1)
     ir2=2*(subregions(isubregion,6) - 1)
     dir=2

     iax=1
     iay=1
     iar=1

     num_material = subregions(isubregion,7)


  end subroutine define_model_regions



  subroutine define_mesh_regions(USE_REGULAR_MESH,isubregion,NER,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,iproc_xi,iproc_eta,&
       nblayers,ner_layer,ndoublings,ner_doublings,&
       iaddx,iaddy,iaddz,ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar)

    implicit none

    include "constants.h"

    logical USE_REGULAR_MESH

    integer isubregion
    integer NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER
    integer iproc_xi,iproc_eta

    integer ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir
    integer iax,iay,iar
    integer nblayers
    !integer num_material
    integer ndoublings

! topology of the elements
    integer iaddx(NGNOD_EIGHT_CORNERS)
    integer iaddy(NGNOD_EIGHT_CORNERS)
    integer iaddz(NGNOD_EIGHT_CORNERS)
    integer ner_layer(nblayers)
    integer ner_doublings(2)

    ! to avoid compiler warnings
    integer idummy
    idummy = ner_layer(1)
    idummy = iproc_xi
    idummy = iproc_eta

! **************

!
!--- case of a regular mesh
!
  if(USE_REGULAR_MESH) then

     call usual_hex_nodes(NGNOD_EIGHT_CORNERS,iaddx,iaddy,iaddz)


     ix1=0
     ix2=2*(NEX_PER_PROC_XI - 1)
     dix=2

     iy1=0
     iy2=2*(NEX_PER_PROC_ETA - 1)
     diy=2

     ir1=0
     ir2=2*(NER - 1)
     dir=2

     iax=1
     iay=1
     iar=1

  else
     if(ndoublings == 1) then

        select case (isubregion)

        case (1)
           ix1=0
           ix2=2*(NEX_PER_PROC_XI - 1)
           dix=4

           iy1=0
           iy2=2*(NEX_PER_PROC_ETA - 1)
           diy=4

           ir1=0
           ir2=2*(ner_doublings(1) -2 -1)
           dir=2

           iax=2
           iay=2
           iar=1

        case (2)

           ix1=0
           ix2=2*(NEX_PER_PROC_XI - 1)
           dix=8

           iy1=0
           iy2=2*(NEX_PER_PROC_ETA - 1)
           diy=8

           ir1=2*(ner_doublings(1) - 2)
           ir2=2*(ner_doublings(1) - 2)
           dir=2

           iax=4
           iay=4
           iar=2

        case (3)

           ix1=0
           ix2=2*(NEX_PER_PROC_XI - 1)
           dix=2

           iy1=0
           iy2=2*(NEX_PER_PROC_ETA - 1)
           diy=2

           ir1=2*(ner_doublings(1))
           ir2=2*(NER - 1)
           dir=2

           iax=1
           iay=1
           iar=1

        case default
           stop 'Wrong number of subregions'

        end select

     else if(ndoublings == 2) then

        select case (isubregion)

        case (1)
           ix1=0
           ix2=2*(NEX_PER_PROC_XI - 1)
           dix=8

           iy1=0
           iy2=2*(NEX_PER_PROC_ETA - 1)
           diy=8

           ir1=0
           ir2=2*(ner_doublings(2) -2 -1)
           dir=2

           iax=4
           iay=4
           iar=1

        case (2)

           ix1=0
           ix2=2*(NEX_PER_PROC_XI - 1)
           dix=16

           iy1=0
           iy2=2*(NEX_PER_PROC_ETA - 1)
           diy=16

           ir1=2*(ner_doublings(2) - 2)
           ir2=2*(ner_doublings(2) - 2)
           dir=2

           iax=8
           iay=8
           iar=2


        case (3)

           ix1=0
           ix2=2*(NEX_PER_PROC_XI - 1)
           dix=4

           iy1=0
           iy2=2*(NEX_PER_PROC_ETA - 1)
           diy=4

           ir1=2*ner_doublings(2)
           ir2=2*(ner_doublings(1) -2 -1)
           dir=2

           iax=2
           iay=2
           iar=1

        case (4)

           ix1=0
           ix2=2*(NEX_PER_PROC_XI - 1)
           dix=8

           iy1=0
           iy2=2*(NEX_PER_PROC_ETA - 1)
           diy=8

           ir1=2*(ner_doublings(1) - 2)
           ir2=2*(ner_doublings(1) - 2)
           dir=2

           iax=4
           iay=4
           iar=2

        case (5)

           ix1=0
           ix2=2*(NEX_PER_PROC_XI - 1)
           dix=2

           iy1=0
           iy2=2*(NEX_PER_PROC_ETA - 1)
           diy=2

           ir1=2*(ner_doublings(1))
           ir2=2*(NER - 1)
           dir=2

           iax=1
           iay=1
           iar=1

        case default
           stop 'Wrong number of subregions'

        end select


     else
        stop 'Wrong number of doublings'

     endif

  endif

  end subroutine define_mesh_regions

