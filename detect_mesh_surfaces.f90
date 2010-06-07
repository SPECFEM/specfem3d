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
  use specfem_par_movie
  use specfem_par_acoustic
  use specfem_par_elastic
  implicit none

  ! for mesh surface
  allocate(ispec_is_surface_external_mesh(NSPEC_AB))
  allocate(iglob_is_surface_external_mesh(NGLOB_AB))

! determines model surface  
  if (.not. RECVS_CAN_BE_BURIED_EXT_MESH .or. &
      EXTERNAL_MESH_MOVIE_SURFACE .or. EXTERNAL_MESH_CREATE_SHAKEMAP) then

    ! returns surface points/elements 
    ! in ispec_is_surface_external_mesh / iglob_is_surface_external_mesh and
    ! number of faces in nfaces_surface_ext_mesh
    call detect_surface(NPROC,NGLOB_AB,NSPEC_AB,ibool,&
                      ispec_is_surface_external_mesh, &
                      iglob_is_surface_external_mesh, &
                      nfaces_surface_ext_mesh, &
                      num_interfaces_ext_mesh, &
                      max_nibool_interfaces_ext_mesh, &
                      nibool_interfaces_ext_mesh, &
                      my_neighbours_ext_mesh, &
                      ibool_interfaces_ext_mesh) 
  endif 

! takes cross-section surfaces instead
  if( (EXTERNAL_MESH_MOVIE_SURFACE .or. EXTERNAL_MESH_CREATE_SHAKEMAP) &
     .and. PLOT_CROSS_SECTIONS ) then
    call detect_surface_cross_section(NPROC,NGLOB_AB,NSPEC_AB,ibool,&
                            ispec_is_surface_external_mesh, &
                            iglob_is_surface_external_mesh, &
                            nfaces_surface_ext_mesh, &
                            num_interfaces_ext_mesh, &
                            max_nibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh, &
                            my_neighbours_ext_mesh, &
                            ibool_interfaces_ext_mesh,&
                            CROSS_SECTION_X,CROSS_SECTION_Y,CROSS_SECTION_Z, &
                            xstore,ystore,zstore,myrank)    
  endif
  
! takes number of faces for top, free surface only
  if( MOVIE_SURFACE .or. CREATE_SHAKEMAP ) then
    nfaces_surface_ext_mesh = num_free_surface_faces
    ! face corner indices
    iorderi(1) = 1
    iorderi(2) = NGLLX
    iorderi(3) = NGLLX
    iorderi(4) = 1
    iorderj(1) = 1
    iorderj(2) = 1
    iorderj(3) = NGLLY
    iorderj(4) = NGLLY    
  endif
  
! handles movies and shakemaps
  if( EXTERNAL_MESH_MOVIE_SURFACE .or. &
     EXTERNAL_MESH_CREATE_SHAKEMAP .or. &
     MOVIE_SURFACE .or. &
     CREATE_SHAKEMAP ) then
    call setup_movie_meshes()
  endif

! stores wavefields for whole volume
  if (MOVIE_VOLUME) then
    ! acoustic
    if( ACOUSTIC_SIMULATION .or. ELASTIC_SIMULATION ) then  
      allocate(velocity_x(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      allocate(velocity_y(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      allocate(velocity_z(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    endif
    ! elastic only
    if( ELASTIC_SIMULATION ) then
      allocate(div(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      allocate(curl_x(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      allocate(curl_y(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      allocate(curl_z(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      div(:,:,:,:) = 0._CUSTOM_REAL
      curl_x(:,:,:,:) = 0._CUSTOM_REAL
      curl_y(:,:,:,:) = 0._CUSTOM_REAL
      curl_z(:,:,:,:) = 0._CUSTOM_REAL
    endif
  endif

! initializes cross-section gif image
  if( PNM_GIF_IMAGE ) then
    call write_PNM_GIF_initialize()
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

  end subroutine detect_mesh_surfaces
  
  
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

