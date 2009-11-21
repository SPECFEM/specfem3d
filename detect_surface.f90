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

  subroutine detect_surface(NPROC,nglob,nspec,ibool,&
                            ispec_is_surface_external_mesh, &
                            iglob_is_surface_external_mesh, &
                            nfaces_surface_external_mesh, &
                            num_interfaces_ext_mesh, &
                            max_nibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh, &
                            my_neighbours_ext_mesh, &
                            ibool_interfaces_ext_mesh)

! detects surface (points/elements) of model based upon valence
!
! returns: ispec_is_surface_external_mesh, iglob_is_surface_external_mesh 
!               and nfaces_surface_external_mesh

  implicit none
  
  include "constants.h"
  
! global indexing  
  integer :: NPROC,nglob,nspec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec):: ibool

! surface  
  logical, dimension(nspec) :: ispec_is_surface_external_mesh
  logical, dimension(nglob) :: iglob_is_surface_external_mesh
  integer :: nfaces_surface_external_mesh

! MPI partitions
  integer :: num_interfaces_ext_mesh
  integer :: max_nibool_interfaces_ext_mesh
  integer,dimension(num_interfaces_ext_mesh):: nibool_interfaces_ext_mesh
  integer,dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh):: ibool_interfaces_ext_mesh
  integer,dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  
!local parameters
  integer, dimension(:), allocatable :: valence_external_mesh
  integer :: ispec,i,j,k,ii,jj,kk,iglob,ier
  
! detecting surface points/elements (based on valence check on NGLL points) for external mesh
  allocate(valence_external_mesh(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocate valence array'

! initialize surface indices
  ispec_is_surface_external_mesh(:) = .false.
  iglob_is_surface_external_mesh(:) = .false.    
  valence_external_mesh(:) = 0

  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          if( iglob < 1 .or. iglob > nglob) then
            print*,'error valence iglob:',iglob,i,j,k,ispec
            stop 'error valence'
          endif
          valence_external_mesh(iglob) = valence_external_mesh(iglob) + 1
        enddo
      enddo
    enddo
  enddo

  ! adds contributions from different partitions to valence_external_mesh
  call assemble_MPI_scalar_i_ext_mesh(NPROC,nglob,valence_external_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)

! determines spectral elements containing surface points
  do ispec = 1, nspec

    ! loops over GLL points not on edges or corners
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          if ( &
           (k == 1 .or. k == NGLLZ) .and. (j /= 1 .and. j /= NGLLY) .and. (i /= 1 .and. i /= NGLLX) .or. &
           (j == 1 .or. j == NGLLY) .and. (k /= 1 .and. k /= NGLLZ) .and. (i /= 1 .and. i /= NGLLX) .or. &
           (i == 1 .or. i == NGLLX) .and. (k /= 1 .and. k /= NGLLZ) .and. (j /= 1 .and. j /= NGLLY) &
           ) then
            iglob = ibool(i,j,k,ispec)
            if (valence_external_mesh(iglob) == 1) then
              ispec_is_surface_external_mesh(ispec) = .true.

              ! sets flags for all gll points on this face
              if (k == 1 .or. k == NGLLZ) then
                do jj = 1, NGLLY
                  do ii = 1, NGLLX
                    iglob_is_surface_external_mesh(ibool(ii,jj,k,ispec)) = .true.
                  enddo
                enddo
              endif
              if (j == 1 .or. j == NGLLY) then
                do kk = 1, NGLLZ
                  do ii = 1, NGLLX
                    iglob_is_surface_external_mesh(ibool(ii,j,kk,ispec)) = .true.
                  enddo
                enddo
              endif
              if (i == 1 .or. i == NGLLX) then
                do kk = 1, NGLLZ
                  do jj = 1, NGLLY
                    iglob_is_surface_external_mesh(ibool(i,jj,kk,ispec)) = .true.
                  enddo
                enddo
              endif
            endif

          endif
        enddo
      enddo
    enddo

  enddo ! nspec

! counts faces for external-mesh movies and shakemaps
  nfaces_surface_external_mesh = 0
  do ispec = 1, nspec
    iglob = ibool(2,2,1,ispec)
    if (iglob_is_surface_external_mesh(iglob)) then
      nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
    endif
    iglob = ibool(2,2,NGLLZ,ispec)
    if (iglob_is_surface_external_mesh(iglob)) then
      nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
    endif
    iglob = ibool(2,1,2,ispec)
    if (iglob_is_surface_external_mesh(iglob)) then
      nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
    endif
    iglob = ibool(2,NGLLY,2,ispec)
    if (iglob_is_surface_external_mesh(iglob)) then
      nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
    endif
    iglob = ibool(1,2,2,ispec)
    if (iglob_is_surface_external_mesh(iglob)) then
      nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
    endif
    iglob = ibool(NGLLX,2,2,ispec)
    if (iglob_is_surface_external_mesh(iglob)) then
      nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
    endif
  enddo 

  end subroutine detect_surface

!
!-------------------------------------------------------------------------------------------------
!

  subroutine detect_surface_cross_section(NPROC,nglob,nspec,ibool,&
                            ispec_is_surface_external_mesh, &
                            iglob_is_surface_external_mesh, &
                            nfaces_surface_external_mesh, &
                            num_interfaces_ext_mesh, &
                            max_nibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh, &
                            my_neighbours_ext_mesh, &
                            ibool_interfaces_ext_mesh,&
                            x_section,y_section,z_section, &
                            xstore,ystore,zstore,myrank)

! instead of surface of model, this returns cross-section surfaces through model 
! at specified x,y,z - coordinates
!
! note: x,y,z coordinates must coincide with the element (outer-)faces, no planes inside element are taken
!         (this is only a quick & dirty cross-section implementation, no sophisticated interpolation of points considered...)
!
! returns: ispec_is_surface_external_mesh, iglob_is_surface_external_mesh 
!               and nfaces_surface_external_mesh

  implicit none
  
  include "constants.h"
  
! global indexing  
  integer :: NPROC,nglob,nspec,myrank
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec):: ibool

! surface  
  logical, dimension(nspec) :: ispec_is_surface_external_mesh
  logical, dimension(nglob) :: iglob_is_surface_external_mesh
  integer :: nfaces_surface_external_mesh

! MPI partitions
  integer :: num_interfaces_ext_mesh
  integer :: max_nibool_interfaces_ext_mesh
  integer,dimension(num_interfaces_ext_mesh):: nibool_interfaces_ext_mesh
  integer,dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh):: ibool_interfaces_ext_mesh
  integer,dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh

! specified x,y,z - coordinates
  real(kind=CUSTOM_REAL):: x_section,y_section,z_section

! mesh global point coordinates
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore
  
!local parameters
  real(kind=CUSTOM_REAL) :: mindist
  integer, dimension(:), allocatable :: valence_external_mesh
  integer :: ispec,i,j,k,ii,jj,kk,iglob,ier,count
  logical :: has_face
  
! detecting surface points/elements (based on valence check on NGLL points) for external mesh
  allocate(valence_external_mesh(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocate valence array'

! initialize surface indices
  ispec_is_surface_external_mesh(:) = .false.
  iglob_is_surface_external_mesh(:) = .false.    
  valence_external_mesh(:) = 0

! an estimation of the minimum distance between global points
  mindist = minval( (xstore(ibool(1,1,1,:)) - xstore(ibool(2,1,1,:)))**2 &
                  + (ystore(ibool(1,1,1,:)) - ystore(ibool(2,1,1,:)))**2 &
                  + (zstore(ibool(1,1,1,:)) - zstore(ibool(2,1,1,:)))**2 )
  mindist = sqrt(mindist)
  
! sets valence to corresponding to process rank  for points on cross-sections
  count = 0
  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)

          ! x cross-section  
          if( abs( xstore(iglob) - x_section ) < 0.5*mindist ) then
            ! sets valence to 1 for points on cross-sections
            valence_external_mesh(iglob) = myrank+1
            count = count + 1
          endif

          ! y cross-section  
          if( abs( ystore(iglob) - y_section ) < 0.5*mindist ) then
            ! sets valence to 1 for points on cross-sections
            valence_external_mesh(iglob) = myrank+1
            count = count + 1
          endif
          
          ! z cross-section  
          if( abs( zstore(iglob) - z_section ) < 0.5*mindist ) then
            ! sets valence to 1 for points on cross-sections
            valence_external_mesh(iglob) = myrank+1
            count = count + 1
          endif
          
        enddo
      enddo
    enddo
  enddo

! adds contributions from different partitions to valence_external_mesh
  call assemble_MPI_scalar_i_ext_mesh(NPROC,nglob,valence_external_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)


! determines spectral elements containing surface points
! (only counts element outer faces, no planes inside element)
  count = 0
  do ispec = 1, nspec

    ! loops over GLL points not on edges or corners, but inside faces
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          if ( ((k == 1 .or. k == NGLLZ) .and. (j == 2 .and. i == 2)) .or. &
              ((j == 1 .or. j == NGLLY) .and. (k == 2 .and. i == 2)) .or. &
              ((i == 1 .or. i == NGLLX) .and. (k == 2 .and. j == 2)) ) then
           
            iglob = ibool(i,j,k,ispec)
           
            ! considers only points in same process or, if point is shared between two processes, 
            ! only with higher process ranks than itself
            if (valence_external_mesh(iglob) == myrank+1 .or. valence_external_mesh(iglob) > 2*(myrank+1) ) then
            
              has_face = .false.
              

              ! sets flags for all gll points on a face and makes sure it's not inside the element
              ! zmin & zmax face
              if ((k == 1 .or. k == NGLLZ) .and. valence_external_mesh(ibool(3,3,k,ispec)) >= 1 ) then
                has_face = .true.
                do jj = 1, NGLLY
                  do ii = 1, NGLLX
                    iglob_is_surface_external_mesh(ibool(ii,jj,k,ispec)) = .true.
                    ! resets valence to count face only once
                    valence_external_mesh(ibool(ii,jj,k,ispec)) = 0
                  enddo
                enddo                
              endif
              
              ! ymin & ymax
              if ((j == 1 .or. j == NGLLY) .and. valence_external_mesh(ibool(3,j,3,ispec)) >= 1) then
                has_face = .true.
                do kk = 1, NGLLZ
                  do ii = 1, NGLLX
                    iglob_is_surface_external_mesh(ibool(ii,j,kk,ispec)) = .true.
                    ! resets valence to count face only once
                    valence_external_mesh(ibool(ii,j,kk,ispec)) = 0 
                  enddo
                enddo
              endif
              
              ! xmin & xmax
              if ((i == 1 .or. i == NGLLX) .and. valence_external_mesh(ibool(i,3,3,ispec)) >= 1) then
                has_face = .true.
                do kk = 1, NGLLZ
                  do jj = 1, NGLLY
                    iglob_is_surface_external_mesh(ibool(i,jj,kk,ispec)) = .true.
                    ! resets valence to count face only once
                    valence_external_mesh(ibool(i,jj,kk,ispec)) = 0
                  enddo
                enddo
              endif


              ! sets flag for element
              if( has_face ) then
                ispec_is_surface_external_mesh(ispec) = .true.
                count = count+1
              endif

            endif            
          endif
        enddo
      enddo
    enddo

  enddo ! nspec

! counts faces for external-mesh movies and shakemaps
  nfaces_surface_external_mesh = 0
  do ispec = 1, nspec
    if( ispec_is_surface_external_mesh(ispec) ) then
      ! zmin face
      if (iglob_is_surface_external_mesh(ibool(2,2,1,ispec))) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
      endif
      ! zmax
      if (iglob_is_surface_external_mesh(ibool(2,2,NGLLZ,ispec))) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
      endif
      ! ymin 
      if (iglob_is_surface_external_mesh(ibool(2,1,2,ispec))) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
      endif
      ! ymax 
      if (iglob_is_surface_external_mesh(ibool(2,NGLLY,2,ispec))) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
      endif
      !xmin 
      if (iglob_is_surface_external_mesh(ibool(1,2,2,ispec))) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
      endif
      !xmax 
      if (iglob_is_surface_external_mesh(ibool(NGLLX,2,2,ispec))) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
      endif
    endif
  enddo 

  end subroutine detect_surface_cross_section


  