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
  integer, dimension(:,:), allocatable :: buffer_send_scalar_i_ext_mesh
  integer, dimension(:,:), allocatable :: buffer_recv_scalar_i_ext_mesh  
  integer, dimension(:), allocatable :: request_send_scalar_ext_mesh
  integer, dimension(:), allocatable :: request_recv_scalar_ext_mesh  
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

  allocate(buffer_send_scalar_i_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  allocate(buffer_recv_scalar_i_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  allocate(request_send_scalar_ext_mesh(num_interfaces_ext_mesh))
  allocate(request_recv_scalar_ext_mesh(num_interfaces_ext_mesh))

  ! adds contributions from different partitions to valence_external_mesh
  call assemble_MPI_scalar_i_ext_mesh(NPROC,nglob,valence_external_mesh, &
                        buffer_send_scalar_i_ext_mesh,buffer_recv_scalar_i_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,my_neighbours_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

  deallocate(buffer_send_scalar_i_ext_mesh)
  deallocate(buffer_recv_scalar_i_ext_mesh)
  deallocate(request_send_scalar_ext_mesh)
  deallocate(request_recv_scalar_ext_mesh)

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

! counts faces for movies and shakemaps
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
  