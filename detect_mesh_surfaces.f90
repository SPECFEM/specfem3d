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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine detect_mesh_surfaces()

  use specfem_par
  implicit none

! detecting surface points/elements (based on valence check on NGLL points) for external mesh


  allocate(ispec_is_surface_external_mesh(NSPEC_AB))
  allocate(iglob_is_surface_external_mesh(NGLOB_AB))
!  allocate(valence_external_mesh(NGLOB_AB))

  if (.not. RECVS_CAN_BE_BURIED_EXT_MESH .or. EXTERNAL_MESH_MOVIE_SURFACE .or. EXTERNAL_MESH_CREATE_SHAKEMAP) then

    ! returns surace points/elements
    call detect_surface(NPROC,NGLOB_AB,NSPEC_AB,ibool,&
                      ispec_is_surface_external_mesh, &
                      iglob_is_surface_external_mesh, &
                      nfaces_surface_external_mesh, &
                      num_interfaces_ext_mesh, &
                      max_nibool_interfaces_ext_mesh, &
                      nibool_interfaces_ext_mesh, &
                      my_neighbours_ext_mesh, &
                      ibool_interfaces_ext_mesh) 

!
!    valence_external_mesh(:) = 0
!    ispec_is_surface_external_mesh(:) = .false.
!    iglob_is_surface_external_mesh(:) = .false.
!    do ispec = 1, NSPEC_AB
!      do k = 1, NGLLZ
!        do j = 1, NGLLY
!          do i = 1, NGLLX
!            iglob = ibool(i,j,k,ispec)
!            valence_external_mesh(iglob) = valence_external_mesh(iglob) + 1
!          enddo
!        enddo
!      enddo
!    enddo
!
!    allocate(buffer_send_scalar_i_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
!    allocate(buffer_recv_scalar_i_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
!
!    ! adds contributions from different partitions to valence_external_mesh
!    call assemble_MPI_scalar_i_ext_mesh(NPROC,NGLOB_AB,valence_external_mesh, &
!         buffer_send_scalar_i_ext_mesh,buffer_recv_scalar_i_ext_mesh, &
!         num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
!         nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,my_neighbours_ext_mesh, &
!         request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
!
!    do ispec = 1, NSPEC_AB
!      do k = 1, NGLLZ
!        do j = 1, NGLLY
!          do i = 1, NGLLX
!            if ( &
!             (k == 1 .or. k == NGLLZ) .and. (j /= 1 .and. j /= NGLLY) .and. (i /= 1 .and. i /= NGLLX) .or. &
!             (j == 1 .or. j == NGLLY) .and. (k /= 1 .and. k /= NGLLZ) .and. (i /= 1 .and. i /= NGLLX) .or. &
!             (i == 1 .or. i == NGLLX) .and. (k /= 1 .and. k /= NGLLZ) .and. (j /= 1 .and. j /= NGLLY) &
!             ) then
!              iglob = ibool(i,j,k,ispec)
!              if (valence_external_mesh(iglob) == 1) then
!                ispec_is_surface_external_mesh(ispec) = .true.
!
!                if (k == 1 .or. k == NGLLZ) then
!                  do jj = 1, NGLLY
!                    do ii = 1, NGLLX
!                      iglob_is_surface_external_mesh(ibool(ii,jj,k,ispec)) = .true.
!                    enddo
!                  enddo
!                endif
!                if (j == 1 .or. j == NGLLY) then
!                  do kk = 1, NGLLZ
!                    do ii = 1, NGLLX
!                      iglob_is_surface_external_mesh(ibool(ii,j,kk,ispec)) = .true.
!                    enddo
!                  enddo
!                endif
!                if (i == 1 .or. i == NGLLX) then
!                  do kk = 1, NGLLZ
!                    do jj = 1, NGLLY
!                      iglob_is_surface_external_mesh(ibool(i,jj,kk,ispec)) = .true.
!                    enddo
!                  enddo
!                endif
!              endif
!
!            endif
!          enddo
!        enddo
!      enddo
!
!    enddo ! nspec
    
  endif 
  
  ! handles movies and shakemaps
  if (EXTERNAL_MESH_MOVIE_SURFACE .or. EXTERNAL_MESH_CREATE_SHAKEMAP) then
    call setup_movie_meshes()
  endif

!!!! NL NL REGOLITH : runs at cines for asteroid simulations. Elements in contact with surface are part of the regolith layer.
!!$  allocate(ispec_is_regolith(NSPEC_AB))
!!$  ispec_is_regolith(:) = .false.
!!$  do ispec = 1, NSPEC_AB
!!$    do k = 1, NGLLZ
!!$      do j = 1, NGLLY
!!$        do i = 1, NGLLX
!!$          iglob = ibool(i,j,k,ispec)
!!$          if (iglob_is_surface_external_mesh(iglob)) then
!!$            ispec_is_regolith(ispec) = .true.
!!$          endif
!!$        enddo
!!$      enddo
!!$    enddo
!!$  enddo
!!$
!!$  do ispec = 1, NSPEC_AB
!!$    if (ispec_is_regolith(ispec)) then
!!$      do k = 1, NGLLZ
!!$        do j = 1, NGLLY
!!$          do i = 1, NGLLX
!!$             kappastore(i,j,k,ispec) = materials_ext_mesh(1,2)* &
!!$                  (materials_ext_mesh(2,2)*materials_ext_mesh(2,2) - &
!!$                  4.d0*materials_ext_mesh(3,2)*materials_ext_mesh(3,2)/3.d0)
!!$             mustore(i,j,k,ispec) = materials_ext_mesh(1,2)*materials_ext_mesh(3,2)*&
!!$                  materials_ext_mesh(3,2)
!!$
!!$          enddo
!!$        enddo
!!$      enddo
!!$    endif
!!$  enddo
!!$
!!$
!!$  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
!!$  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
!!$  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
!!$
!!$  rmass(:) = 0._CUSTOM_REAL
!!$
!!$  do ispec=1,NSPEC_AB
!!$  do k=1,NGLLZ
!!$    do j=1,NGLLY
!!$      do i=1,NGLLX
!!$        weight=wxgll(i)*wygll(j)*wzgll(k)
!!$        iglob=ibool(i,j,k,ispec)
!!$
!!$        jacobianl=jacobian(i,j,k,ispec)
!!$
!!$! distinguish between single and double precision for reals
!!$        if (.not. ispec_is_regolith(ispec)) then
!!$        if(CUSTOM_REAL == SIZE_REAL) then
!!$          rmass(iglob) = rmass(iglob) + &
!!$               sngl(dble(materials_ext_mesh(1,1)) * dble(jacobianl) * weight)
!!$        else
!!$          rmass(iglob) = rmass(iglob) + materials_ext_mesh(1,1) * jacobianl * weight
!!$        endif
!!$        else
!!$        if(CUSTOM_REAL == SIZE_REAL) then
!!$          rmass(iglob) = rmass(iglob) + &
!!$               sngl(dble(materials_ext_mesh(1,2)) * dble(jacobianl) * weight)
!!$        else
!!$          rmass(iglob) = rmass(iglob) + materials_ext_mesh(1,2) * jacobianl * weight
!!$        endif
!!$        endif
!!$
!!$      enddo
!!$    enddo
!!$  enddo
!!$  enddo


!!!! NL NL REGOLITH

!!!!!!!!!! DK DK   endif

  end subroutine
  
  
!!!! NL NL REGOLITH
!!$  double precision function materials_ext_mesh(i,j)
!!$
!!$    implicit none
!!$
!!$    integer :: i,j
!!$
!!$    select case (j)
!!$      case (1)
!!$        select case (i)
!!$          case (1)
!!$            materials_ext_mesh = 2700.d0
!!$          case (2)
!!$            materials_ext_mesh = 3000.d0
!!$          case (3)
!!$            materials_ext_mesh = 1732.051d0
!!$          case default
!!$            call stop_all()
!!$          end select
!!$      case (2)
!!$        select case (i)
!!$          case (1)
!!$            materials_ext_mesh = 2000.d0
!!$          case (2)
!!$            materials_ext_mesh = 900.d0
!!$          case (3)
!!$            materials_ext_mesh = 500.d0
!!$          case default
!!$            call stop_all()
!!$          end select
!!$      case default
!!$        call stop_all()
!!$    end select
!!$
!!$  end function materials_ext_mesh
!!!! NL NL REGOLITH

